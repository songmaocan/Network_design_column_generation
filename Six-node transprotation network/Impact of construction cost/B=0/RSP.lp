\ Model RSP
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  3 x_0_1 + 3 x_0_2 + 3 x_1_0 + 3 x_1_3 + 3 x_2_0 + 3 x_2_3 + 3 x_3_1
   + 3 x_3_2 + 2 x_0_4 + 3.82842712474619 x_1_4 + x_4_5 + x_5_2 + 2 x_5_3
   + y
Subject To
 Node_0: x_0_1 + x_0_2 - x_1_0 - x_2_0 + x_0_4 = 0
 Node_1: - x_0_1 + x_1_0 + x_1_3 - x_3_1 + x_1_4 = 1
 Node_2: - x_0_2 + x_2_0 + x_2_3 - x_3_2 - x_5_2 = -1
 Node_3: - x_1_3 - x_2_3 + x_3_1 + x_3_2 - x_5_3 = 0
 Node_4: - x_0_4 - x_1_4 + x_4_5 = 0
 Node_5: - x_4_5 + x_5_2 + x_5_3 = 0
 variance_limit: 4 x_0_1 + 4 x_0_2 + 4 x_1_0 + 4 x_1_3 + 4 x_2_0 + 4 x_2_3
   + 4 x_3_1 + 4 x_3_2 + x_0_4 + 4 x_1_4 + x_4_5 + 4 x_5_2 + x_5_3 + [
   - y ^2 ] <= 0
Bounds
Binaries
 x_0_1 x_0_2 x_1_0 x_1_3 x_2_0 x_2_3 x_3_1 x_3_2 x_0_4 x_1_4 x_4_5 x_5_2
 x_5_3
End
