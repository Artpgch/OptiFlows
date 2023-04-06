# OptiFlows
Stream optimization algorithm in data transmission networks
This software package consists of files OptiFlows.m and TargetFunction.m.
File OptiFlows.m contains main algorithm code, file TargetFunction.m contains code, which forms an objective function. Initial data and the result of calculations are represented in xlsx-spreadsheets for convenient visualization. Source file NetMatrix.xlsx contains a spreadsheet with a structure of the network, source file FlowMatrix.xlsx contains a spreadsheet with structure of direction and intensity of the flows.
According to the results of calculation, software makes files OptiFlows.xlsx and LoadCoeffs.xlsx.
Result file OptiFlows.xlsx contain a spreadsheets with optimal subflows allocation data.
Result file LoadCoeffs.xlsx contain a spreadsheets with communication channels loading coefficients.
If optimization impossible, the software will not create file OptiFlows.xlsx, but will create file LoadCoeffs.xlsx, which contains data representing in how much times it is necessary to increase certain communication channels for successful optimization.

Files NetMatrix.xlsx and FlowMatrix.xlsx in this repository are examples of the source data.
