[XYGoal, XYObstacles, RadiusObstacles, FieldSize] = LoadRoboCupField();

[Path PathAPF] = SolvePathPlan_LS([20 50] , XYGoal, XYObstacles, RadiusObstacles, FieldSize);
figure; plot(PathAPF);

[Path PathAPF] = SolvePathPlan_LS([100 140] , XYGoal, XYObstacles, RadiusObstacles, FieldSize);
figure; plot(PathAPF);

[Path PathAPF] = SolvePathPlan_LS([120 10] , XYGoal, XYObstacles, RadiusObstacles, FieldSize);
figure; plot(PathAPF);







