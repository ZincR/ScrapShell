function cat = model_catalog()
%MODEL_CATALOG The example shells shown in the viewer (all from the official PH-CPF repo).
ids    = {'metro', 'intersection', 'six_s0', 'vase', '3holes_r', 'car_r', 'sphere2con_s0', 'c_s0'};
labels = {'Metro', 'Intersection', 'Six', 'Vase', 'Three holes', 'Car', 'Sphere with cones', 'C'};
groups = [repmat({'Open shells'}, 1, 4), repmat({'Closed models'}, 1, 4)];
cat = struct('id', ids, 'label', labels, 'group', groups);
end
