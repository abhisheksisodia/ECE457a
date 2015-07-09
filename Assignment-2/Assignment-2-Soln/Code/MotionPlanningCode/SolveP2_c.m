[XYGoal, XYObstacles, RadiusObstacles, FieldSize] = LoadRoboCupField();

SAParams.ITemp = 500;
SAParams.FTemp = 10;
SAParams.Alpha = 0.85;
SAParams.NumItsPerTemp = 10;
        
[Path PathAPF] = SolvePathPlan_LS_SA([20 50] , XYGoal, XYObstacles, RadiusObstacles, FieldSize, SAParams);
size(Path, 1)
figure; plot(PathAPF);

[Path PathAPF] = SolvePathPlan_LS_SA([100 140] , XYGoal, XYObstacles, RadiusObstacles, FieldSize, SAParams);
size(Path, 1)
figure; plot(PathAPF);

[Path PathAPF] = SolvePathPlan_LS_SA([120 10] , XYGoal, XYObstacles, RadiusObstacles, FieldSize, SAParams);
size(Path, 1)
figure; plot(PathAPF);


