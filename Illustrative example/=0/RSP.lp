\ Model RSP
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x_0_1 + 2 x_0_2 + x_1_0 + 2 x_1_3 + 2 x_2_0 + 2 x_2_3 + 2 x_3_1 + x_3_2
   + 2 x_0_3 + 2 x_1_2
Subject To
 Node_0: x_0_1 + x_0_2 - x_1_0 - x_2_0 + x_0_3 = 0
 Node_1: - x_0_1 + x_1_0 + x_1_3 - x_3_1 + x_1_2 = 1
 Node_2: - x_0_2 + x_2_0 + x_2_3 - x_3_2 - x_1_2 = -1
 Node_3: - x_1_3 - x_2_3 + x_3_1 + x_3_2 - x_0_3 = 0
 variance_limit: x_0_1 + 3 x_0_2 + x_1_0 + 3 x_1_3 + x_2_0 + x_2_3 + x_3_1
   + x_3_2 + x_0_3 + 3 x_1_2 + [ - y ^2 ] <= 0
Bounds
Binaries
 x_0_1 x_0_2 x_1_0 x_1_3 x_2_0 x_2_3 x_3_1 x_3_2 x_0_3 x_1_2
End
