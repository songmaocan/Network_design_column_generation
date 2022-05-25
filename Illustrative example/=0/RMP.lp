\ Model Relaxed master problem
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  3 x_0_0 + 2 x_0_1 + 3 x_0_2 + 3 x_0_3 + 3 x_1_0 + x_1_1 + 3 x_1_2 + x_1_3
Subject To
 flow_balance_for_OD-0: x_0_0 + x_0_1 + x_0_2 + x_0_3 = 1
 flow_balance_for_OD-1: x_1_0 + x_1_1 + x_1_2 + x_1_3 = 1
 flow_balance_for_construction: y_0 + y_1 + y_2 + y_3 = 1
 Coupling_0_0: x_0_1 - y_3 <= 0
 Coupling_1_0: - y_2 <= 0
 Coupling_0_1: - y_3 <= 0
 Coupling_1_1: x_1_1 + x_1_3 - y_2 <= 0
Bounds
 x_0_0 <= 1
 x_0_1 <= 1
 x_0_2 <= 1
 x_0_3 <= 1
 x_1_0 <= 1
 x_1_1 <= 1
 x_1_2 <= 1
 x_1_3 <= 1
 y_0 <= 1
 y_1 <= 1
 y_2 <= 1
 y_3 <= 1
End
