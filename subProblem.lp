\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  0.4186046511627907 x0_5_0 + 0.4186046511627907 x0_5_1
   + 0.4186046511627907 x0_5_2 + 0.2790697674418605 x2_4_0
   + 0.2790697674418605 x2_4_1 + 0.2790697674418605 x2_4_2
   + 0.3023255813953488 x2_6_0 + 0.3023255813953488 x2_6_1
   + 0.3023255813953488 x2_6_2 - 669.7674418604652 Constant
Subject To
 Supply_cons0_0: x0_0_0 + x0_1_0 + x0_2_0 + x0_3_0 + x0_4_0 + x0_5_0
   + x0_6_0 = 400
 Supply_cons0_1: x0_0_1 + x0_1_1 + x0_2_1 + x0_3_1 + x0_4_1 + x0_5_1
   + x0_6_1 = 800
 Supply_cons0_2: x0_0_2 + x0_1_2 + x0_2_2 + x0_3_2 + x0_4_2 + x0_5_2
   + x0_6_2 = 200
 Supply_cons1_0: x1_0_0 + x1_1_0 + x1_2_0 + x1_3_0 + x1_4_0 + x1_5_0
   + x1_6_0 = 700
 Supply_cons1_1: x1_0_1 + x1_1_1 + x1_2_1 + x1_3_1 + x1_4_1 + x1_5_1
   + x1_6_1 = 1600
 Supply_cons1_2: x1_0_2 + x1_1_2 + x1_2_2 + x1_3_2 + x1_4_2 + x1_5_2
   + x1_6_2 = 300
 Supply_cons2_0: x2_0_0 + x2_1_0 + x2_2_0 + x2_3_0 + x2_4_0 + x2_5_0
   + x2_6_0 = 800
 Supply_cons2_1: x2_0_1 + x2_1_1 + x2_2_1 + x2_3_1 + x2_4_1 + x2_5_1
   + x2_6_1 = 1800
 Supply_cons2_2: x2_0_2 + x2_1_2 + x2_2_2 + x2_3_2 + x2_4_2 + x2_5_2
   + x2_6_2 = 300
 Demand_cons0_0: x0_0_0 + x1_0_0 + x2_0_0 = 300
 Demand_cons0_1: x0_0_1 + x1_0_1 + x2_0_1 = 500
 Demand_cons0_2: x0_0_2 + x1_0_2 + x2_0_2 = 100
 Demand_cons1_0: x0_1_0 + x1_1_0 + x2_1_0 = 300
 Demand_cons1_1: x0_1_1 + x1_1_1 + x2_1_1 = 750
 Demand_cons1_2: x0_1_2 + x1_1_2 + x2_1_2 = 100
 Demand_cons2_0: x0_2_0 + x1_2_0 + x2_2_0 = 100
 Demand_cons2_1: x0_2_1 + x1_2_1 + x2_2_1 = 400
 Demand_cons2_2: x0_2_2 + x1_2_2 + x2_2_2 = 0
 Demand_cons3_0: x0_3_0 + x1_3_0 + x2_3_0 = 75
 Demand_cons3_1: x0_3_1 + x1_3_1 + x2_3_1 = 250
 Demand_cons3_2: x0_3_2 + x1_3_2 + x2_3_2 = 50
 Demand_cons4_0: x0_4_0 + x1_4_0 + x2_4_0 = 650
 Demand_cons4_1: x0_4_1 + x1_4_1 + x2_4_1 = 950
 Demand_cons4_2: x0_4_2 + x1_4_2 + x2_4_2 = 200
 Demand_cons5_0: x0_5_0 + x1_5_0 + x2_5_0 = 225
 Demand_cons5_1: x0_5_1 + x1_5_1 + x2_5_1 = 850
 Demand_cons5_2: x0_5_2 + x1_5_2 + x2_5_2 = 100
 Demand_cons6_0: x0_6_0 + x1_6_0 + x2_6_0 = 250
 Demand_cons6_1: x0_6_1 + x1_6_1 + x2_6_1 = 500
 Demand_cons6_2: x0_6_2 + x1_6_2 + x2_6_2 = 250
Bounds
 Constant = 1
End
