// phcpf_mesher: replaces PH-CPF's Windows-only MeshCutter / Parametrization / Mesher executables,
// keeping the original division of work (see CPF.m, uv2mesh / cut_and_comb / param_and_mesh):
//
//   setup  MeshCutter + Parametrization: read the matching and singularities computed in MATLAB
//          (CPF.compute_matching_and_sing), cut and comb, solve the seamless Poisson problem without
//          integers, and export the system that the official BarrierPoissonSolve.m expects
//          (A, b, C, G, x0, integerIndices, x2CornerMat, cut mesh V/F, face normals FN, N).
//   mesh   Mesher: take the rounded solution x from BarrierPoissonSolve and extract the isoline mesh.
//   auto   Directional's own integer rounding, without the barrier (smoke tests only).
//
// Usage: phcpf_mesher <setup|mesh|auto> mesh.off field.rawfield out_prefix
//            [--matching file] [--sing file] [--x file] [--lr L]
//   --matching  lines "v0 v1 fi fj m" (0-based): CPF matching m of the inner edge (v0, v1), mapping
//               vector k of face fi to vector k+m of face fj. Without it, curl matching is used.
//   --sing      lines "v index" (0-based): CPF singularities. Without it, the matching's holonomy.
//   --x         (mesh mode) the rounded variable vector, one value per line.
//   --lr        length ratio; default 1/bbd. Current Directional divides the field by its mean norm and
//               the parameter length by the same factor (integrate.h), so 1/bbd keeps the raw field's
//               own magnitude, which carries the CPF sizing. (The official code passed 1/(navg*bbd) to
//               an older Directional that only normalized the field.)
//
// Outputs (prefix_*): setup -> A.txt b.txt C.txt G.txt x2c.txt (sparse triplets), x0.txt int.txt
// singidx.txt cutV.txt cutF.txt FN.txt meta.json; mesh/auto -> tri.off corner.txt sing.txt stats.json.
#include <chrono>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <directional/TriMesh.h>
#include <directional/PCFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/readOFF.h>
#include <directional/read_raw_field.h>
#include <directional/effort_to_indices.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/branched_gradient.h>
#include <directional/setup_mesher.h>
#include <directional/mesher.h>

using namespace Eigen;
using Clock = std::chrono::steady_clock;
static double secs(Clock::time_point a, Clock::time_point b) { return std::chrono::duration<double>(b - a).count(); }

// ---------------------------------------------------------------------------------------------
// Matching and singularities
// ---------------------------------------------------------------------------------------------

// The official pipeline passes MATLAB's matching to the cutter. Directional's convention:
// matching(e) = j means vector k of EF(e,0) matches vector (k+j) % N of EF(e,1).
static bool import_matching(directional::CartesianField& field, const directional::TriMesh& mesh,
                            const std::string& path) {
  std::ifstream in(path);
  if (!in) return false;
  const int N = field.N;
  std::map<std::pair<int, int>, std::array<int, 3>> byEdge;
  int v0, v1, fi, fj, m;
  while (in >> v0 >> v1 >> fi >> fj >> m) byEdge[{std::min(v0, v1), std::max(v0, v1)}] = {fi, fj, m};
  field.matching.resize(mesh.EF.rows());
  field.matching.setConstant(-1);
  int missing = 0;
  for (int e = 0; e < mesh.EF.rows(); e++) {
    if (mesh.EF(e, 0) == -1 || mesh.EF(e, 1) == -1) continue;
    auto it = byEdge.find({std::min(mesh.EV(e, 0), mesh.EV(e, 1)), std::max(mesh.EV(e, 0), mesh.EV(e, 1))});
    if (it == byEdge.end()) { missing++; continue; }
    auto [a, b, mm] = it->second;
    if (mesh.EF(e, 0) == a && mesh.EF(e, 1) == b) field.matching(e) = ((mm % N) + N) % N;
    else if (mesh.EF(e, 0) == b && mesh.EF(e, 1) == a) field.matching(e) = ((-mm % N) + N) % N;
    else { missing++; }
  }
  if (missing) std::cerr << "phcpf_mesher: " << missing << " inner edges without an imported matching\n";
  field.effort = VectorXd::Zero(mesh.EF.rows());
  return missing == 0;
}

// Fallback when no MATLAB matching is given: per edge, the rotation under which the gradients agree
// along the edge. (directional::curl_matching computes the same matching but indexes the effort with
// `2 * (m + j + N) % N`, i.e. (2*(m+j+N)) % N, so its singularities are wrong; we only need the
// matching and take singularities from the holonomy below.)
static void curl_matching(directional::CartesianField& field, const directional::TriMesh& mesh) {
  const int N = field.N;
  field.matching.resize(mesh.EF.rows());
  field.matching.setConstant(-1);
  field.effort = VectorXd::Zero(mesh.EF.rows());
  for (int i = 0; i < mesh.EF.rows(); i++) {
    if (mesh.EF(i, 0) == -1 || mesh.EF(i, 1) == -1) continue;
    const RowVector3d ev = (mesh.V.row(mesh.EV(i, 1)) - mesh.V.row(mesh.EV(i, 0))).normalized();
    int best = 0; double minCurl = 1e300;
    for (int j = 0; j < N; j++) {
      double c = 0;
      for (int k = 0; k < N; k++) {
        RowVector3d d = field.extField.block(mesh.EF(i, 1), 3 * ((j + k) % N), 1, 3) -
                        field.extField.block(mesh.EF(i, 0), 3 * k, 1, 3);
        c += std::pow(ev.dot(d), 2.0);
      }
      if (c < minCurl) { minCurl = c; best = j; }
    }
    field.matching(i) = best;
  }
}

// Holonomy of the matching around each interior vertex, mod N, in (-N/2, N/2].
static VectorXi holonomy(const directional::CartesianField& field, const directional::TriMesh& mesh) {
  const int N = field.N;
  VectorXi hol = VectorXi::Zero(mesh.V.rows());
  for (int e = 0; e < mesh.EF.rows(); e++) {
    if (mesh.EF(e, 0) == -1 || mesh.EF(e, 1) == -1) continue;
    hol(mesh.EV(e, 0)) -= field.matching(e);
    hol(mesh.EV(e, 1)) += field.matching(e);
  }
  for (int v = 0; v < hol.size(); v++) {
    int h = ((hol(v) % N) + N) % N;
    hol(v) = (mesh.isBoundaryVertex(v) ? 0 : (h > N / 2 ? h - N : h));
  }
  return hol;
}

static void set_singularities(directional::CartesianField& field, const std::vector<int>& v, const std::vector<int>& idx) {
  VectorXi c = Map<const VectorXi>(v.data(), v.size());
  VectorXi s = Map<const VectorXi>(idx.data(), idx.size());
  field.set_singularities(c, s);
}

// ---------------------------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------------------------
static void write_triplets(const std::string& path, const SparseMatrix<double>& M) {
  std::ofstream f(path); f.precision(17);
  f << M.rows() << " " << M.cols() << " " << M.nonZeros() << "\n";
  for (int k = 0; k < M.outerSize(); ++k)
    for (SparseMatrix<double>::InnerIterator it(M, k); it; ++it) f << it.row() << " " << it.col() << " " << it.value() << "\n";
}
template <class Derived> static void write_dense(const std::string& path, const MatrixBase<Derived>& M) {
  std::ofstream f(path); f.precision(17);
  for (int i = 0; i < M.rows(); i++) { for (int j = 0; j < M.cols(); j++) f << (j ? " " : "") << M(i, j); f << "\n"; }
}
static VectorXd read_vector(const std::string& path) {
  std::ifstream f(path); std::vector<double> v; double x;
  while (f >> x) v.push_back(x);
  return Map<VectorXd>(v.data(), v.size());
}

// Rank-reduced constraint matrix, exactly as in directional::integrate.
static SparseMatrix<double> reduce_rows(const SparseMatrix<double>& C) {
  if (C.rows() == 0) return C;
  SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> qr;
  qr.compute(C.transpose());
  const int rank = qr.rank();
  VectorXi P = qr.colsPermutation().indices();
  std::vector<int> pos(C.rows(), -1);
  for (int j = 0; j < rank; j++) pos[P(j)] = j;
  std::vector<Triplet<double>> T;
  for (int k = 0; k < C.outerSize(); ++k)
    for (SparseMatrix<double>::InnerIterator it(C, k); it; ++it)
      if (pos[it.row()] >= 0) T.emplace_back(pos[it.row()], it.col(), it.value());
  SparseMatrix<double> R(rank, C.cols());
  R.setFromTriplets(T.begin(), T.end());
  return R;
}

int main(int argc, char** argv) {
  if (argc < 5) {
    std::cerr << "usage: phcpf_mesher <setup|mesh|auto> mesh.off field.rawfield out_prefix "
                 "[--matching f] [--sing f] [--x f] [--lr L]\n";
    return 2;
  }
  const std::string mode = argv[1], meshPath = argv[2], fieldPath = argv[3], out = argv[4];
  std::string matchingPath, singPath, xPath;
  double lr = -1.0;
  for (int i = 5; i < argc; i++) {
    if (!std::strcmp(argv[i], "--matching") && i + 1 < argc) matchingPath = argv[++i];
    else if (!std::strcmp(argv[i], "--sing") && i + 1 < argc) singPath = argv[++i];
    else if (!std::strcmp(argv[i], "--x") && i + 1 < argc) xPath = argv[++i];
    else if (!std::strcmp(argv[i], "--lr") && i + 1 < argc) lr = std::atof(argv[++i]);
    else { std::cerr << "unknown argument " << argv[i] << "\n"; return 2; }
  }
  if (mode != "setup" && mode != "mesh" && mode != "auto") { std::cerr << "unknown mode " << mode << "\n"; return 2; }
  if (mode == "mesh" && xPath.empty()) { std::cerr << "mesh mode needs --x\n"; return 2; }

  auto t0 = Clock::now();
  directional::TriMesh meshWhole, meshCut;
  directional::PCFaceTangentBundle ftb;
  directional::CartesianField rawField, combedField;
  int N = 0;
  if (!directional::readOFF(meshPath, meshWhole)) { std::cerr << "cannot read " << meshPath << "\n"; return 1; }
  ftb.init(meshWhole);
  if (!directional::read_raw_field(fieldPath, ftb, N, rawField)) { std::cerr << "cannot read " << fieldPath << "\n"; return 1; }
  if (N != 6) { std::cerr << "expected N=6, got " << N << "\n"; return 1; }
  const double bbd = (meshWhole.V.colwise().maxCoeff() - meshWhole.V.colwise().minCoeff()).norm();
  if (lr <= 0) lr = 1.0 / bbd;

  // --- matching and singularities (MATLAB's, as the official cutter received them)
  if (matchingPath.empty()) curl_matching(rawField, meshWhole);
  else if (!import_matching(rawField, meshWhole, matchingPath)) { std::cerr << "bad matching file\n"; return 1; }
  VectorXi hol = holonomy(rawField, meshWhole);
  std::vector<int> sv, si;
  if (!singPath.empty()) {
    std::ifstream in(singPath); int v, s;
    while (in >> v >> s) if (!meshWhole.isBoundaryVertex(v)) { sv.push_back(v); si.push_back(s); }
    int disagree = 0;
    std::vector<bool> isS(meshWhole.V.rows(), false);
    for (size_t k = 0; k < sv.size(); k++) { isS[sv[k]] = true; if (((si[k] - hol(sv[k])) % N) != 0) disagree++; }
    for (int v2 = 0; v2 < hol.size(); v2++) if (hol(v2) != 0 && !isS[v2]) disagree++;
    if (disagree) std::cerr << "phcpf_mesher: warning, " << disagree << " vertices where the singularities and the matching's holonomy disagree\n";
  } else {
    for (int v = 0; v < hol.size(); v++) if (hol(v) != 0) { sv.push_back(v); si.push_back(hol(v)); }
  }
  set_singularities(rawField, sv, si);
  int sumIdx = 0; for (int s : si) sumIdx += s;

  directional::IntegrationData intData(N);
  intData.set_triangular_symmetry(N);
  directional::setup_integration(rawField, intData, meshCut, combedField);
  intData.verbose = false;
  intData.lengthRatio = lr;
  intData.roundSeams = false;         // auto mode: round singularities, then the remaining seams
  auto t1 = Clock::now();
  std::cerr << "phcpf_mesher[" << mode << "]: " << (matchingPath.empty() ? "curl" : "MATLAB") << " matching, "
            << sv.size() << " singularities (index sum " << sumIdx << "/" << N << "), lengthRatio " << lr << "\n";

  // P maps the compressed variables to the N functions at each cut-mesh vertex.
  SparseMatrix<double> P = intData.vertexTrans2CutMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
  VectorXi integerIndices(intData.integerVars.size() * intData.n);
  for (int i = 0; i < intData.integerVars.size(); i++)
    for (int j = 0; j < intData.n; j++) integerIndices(intData.n * i + j) = intData.n * intData.integerVars(i) + j;

  if (mode == "setup") {
    // Continuous seamless solution (no integers), as the official Parametrization_bin without -I.
    directional::IntegrationData cont = intData;
    cont.integralSeamless = false;
    MatrixXd NF, NCF;
    directional::integrate(combedField, cont, meshCut, NF, NCF);
    const VectorXd x0 = cont.nVertexFunction;

    // Poisson system, assembled exactly as in directional::integrate (edge weights are all 1).
    MatrixXd raw = combedField.extField;
    double avg = 0;
    for (int i = 0; i < meshWhole.F.rows(); i++) for (int j = 0; j < N; j++) avg += raw.block(i, 3 * j, 1, 3).norm();
    avg /= double(N * meshWhole.F.rows());
    raw.array() /= avg;
    const double paramLength = bbd * lr / avg;
    std::vector<Triplet<double>> d0T;
    VectorXd gamma(3 * N * meshWhole.F.rows());
    for (int i = 0; i < meshCut.F.rows(); i++)
      for (int j = 0; j < 3; j++)
        for (int k = 0; k < N; k++) {
          d0T.emplace_back(3 * N * i + N * j + k, N * meshCut.F(i, j) + k, -1.0);
          d0T.emplace_back(3 * N * i + N * j + k, N * meshCut.F(i, (j + 1) % 3) + k, 1.0);
          Vector3d ev = (meshCut.V.row(meshCut.F(i, (j + 1) % 3)) - meshCut.V.row(meshCut.F(i, j))).transpose();
          gamma(3 * N * i + N * j + k) = (raw.block(i, 3 * k, 1, 3) * ev)(0, 0) / paramLength;
        }
    SparseMatrix<double> d0(3 * N * meshWhole.F.rows(), N * meshCut.V.rows());
    d0.setFromTriplets(d0T.begin(), d0T.end());
    SparseMatrix<double> A = d0 * P;
    SparseMatrix<double> C = reduce_rows(intData.constraintMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat);
    SparseMatrix<double> G;
    directional::branched_gradient(meshCut, N, G);
    SparseMatrix<double> Gd = G * P;

    write_triplets(out + "_A.txt", A);
    write_dense(out + "_b.txt", gamma);
    write_triplets(out + "_C.txt", C);
    write_triplets(out + "_G.txt", Gd);
    write_triplets(out + "_x2c.txt", P);
    write_dense(out + "_x0.txt", x0);
    write_dense(out + "_int.txt", integerIndices);
    write_dense(out + "_singidx.txt", intData.singularIndices);
    write_dense(out + "_cutV.txt", meshCut.V);
    write_dense(out + "_cutF.txt", meshCut.F);
    write_dense(out + "_FN.txt", meshCut.faceNormals);
    {
      std::ofstream f(out + "_sing.txt");
      for (size_t k = 0; k < sv.size(); k++) f << sv[k] << " " << si[k] << "\n";
    }
    std::ofstream(out + "_meta.json") << "{\"N\":" << N << ",\"n\":" << intData.n << ",\"numVars\":" << P.cols()
                                      << ",\"nInteger\":" << integerIndices.size() << ",\"nSing\":" << sv.size()
                                      << ",\"lengthRatio\":" << lr << ",\"bbd\":" << bbd
                                      << ",\"secSetup\":" << secs(t0, Clock::now()) << "}\n";
    std::cerr << "phcpf_mesher[setup]: " << P.cols() << " variables, " << integerIndices.size()
              << " integer, " << C.rows() << " constraints\n";
    return 0;
  }

  // --- mesh / auto: obtain the (rounded) compressed solution
  MatrixXd NFunction, NCornerFunction;
  if (mode == "auto") {
    intData.integralSeamless = true;
    directional::integrate(combedField, intData, meshCut, NFunction, NCornerFunction);
  } else {
    VectorXd x = read_vector(xPath);
    if (x.size() != P.cols()) { std::cerr << "--x has " << x.size() << " values, expected " << P.cols() << "\n"; return 1; }
    // Snap the rounded variables exactly (BarrierPoissonSolve leaves them within 1e-6); the exact mesher
    // needs singular vertices on lattice points and integer seam translations to match seams.
    for (int i = 0; i < integerIndices.size(); i++) x(integerIndices(i)) = std::round(x(integerIndices(i)));
    for (int i = 0; i < intData.singularIndices.size(); i++)
      x(intData.singularIndices(i)) = std::round(x(intData.singularIndices(i)));
    intData.nVertexFunction = x;
    VectorXd nf = P * x;
    NFunction.resize(meshCut.V.rows(), N);
    for (int i = 0; i < NFunction.rows(); i++) NFunction.row(i) = nf.segment(N * i, N).transpose();
    NCornerFunction.resize(meshWhole.F.rows(), 3 * N);
    for (int i = 0; i < meshWhole.F.rows(); i++)
      for (int j = 0; j < 3; j++) NCornerFunction.block(i, N * j, 1, N) = NFunction.row(meshCut.F(i, j));
  }
  auto t2 = Clock::now();
  {
    std::ofstream f(out + "_corner.txt"); f.precision(10);
    f << NCornerFunction.rows() << " " << N << "\n";
    for (int i = 0; i < NCornerFunction.rows(); i++) {
      for (int j = 0; j < NCornerFunction.cols(); j++) f << (j ? " " : "") << NCornerFunction(i, j);
      f << "\n";
    }
  }
  {
    std::ofstream f(out + "_sing.txt");
    for (size_t k = 0; k < sv.size(); k++) f << sv[k] << " " << si[k] << "\n";
  }

  directional::MesherData mData;
  directional::setup_mesher(meshCut, intData, mData);
  mData.verbose = false;
  MatrixXd V; VectorXi D; MatrixXi F;
  directional::mesher(meshWhole, mData, V, D, F);
  auto t3 = Clock::now();
  {
    std::ofstream f(out + "_tri.off"); f.precision(17);
    f << "OFF\n" << V.rows() << " " << F.rows() << " 0\n";
    for (int i = 0; i < V.rows(); i++) f << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
    for (int i = 0; i < F.rows(); i++) { f << D(i); for (int j = 0; j < D(i); j++) f << " " << F(i, j); f << "\n"; }
  }
  std::map<int, int> hist;
  for (int i = 0; i < D.size(); i++) hist[D(i)]++;
  {
    std::ofstream f(out + "_stats.json");
    f << "{\"nV\":" << V.rows() << ",\"nF\":" << F.rows() << ",\"degreeHist\":{";
    bool first = true;
    for (auto& [d, c] : hist) { f << (first ? "" : ",") << "\"" << d << "\":" << c; first = false; }
    f << "},\"nSing\":" << sv.size() << ",\"lengthRatio\":" << lr << ",\"bbd\":" << bbd
      << ",\"secSetup\":" << secs(t0, t1) << ",\"secSolve\":" << secs(t1, t2) << ",\"secMesh\":" << secs(t2, t3) << "}\n";
  }
  std::cout << "phcpf_mesher[" << mode << "]: " << F.rows() << " polygons, " << V.rows() << " vertices, "
            << sv.size() << " singularities, degrees:";
  for (auto& [d, c] : hist) std::cout << " " << d << ":" << c;
  std::cout << "  (" << secs(t0, t3) << " s)\n";
  return 0;
}
