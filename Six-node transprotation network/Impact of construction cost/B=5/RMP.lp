\ Model Relaxed master problem
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  6.732050807568877 x_0_0 + 6 x_1_0
Subject To
 flow_balance_for_OD-0: x_0_0 = 1
 flow_balance_for_OD-1: x_1_0 = 1
 flow_balance_for_construction: y_0 = 1
 Coupling_0_0: x_0_0 - y_0 <= 0
 Coupling_1_0: - y_0 <= 0
 Coupling_2_0: x_0_0 - y_0 <= 0
 Coupling_3_0: - y_0 <= 0
 Coupling_4_0: x_0_0 - y_0 <= 0
 Coupling_0_1: - y_0 <= 0
 Coupling_1_1: x_1_0 - y_0 <= 0
 Coupling_2_1: x_1_0 - y_0 <= 0
 Coupling_3_1: x_1_0 - y_0 <= 0
 Coupling_4_1: - y_0 <= 0
Bounds
 x_0_0 <= 1
 x_1_0 <= 1
 y_0 <= 1
End
