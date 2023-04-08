\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  Artificial
Subject To
 Capacity_cons0_0: - Artificial + 225 lambda3 <= 625
 Capacity_cons0_1: - Artificial + 650 lambda2 + 150 lambda3 <= 625
 Capacity_cons0_2: - Artificial <= 625
 Capacity_cons0_3: - Artificial <= 625
 Capacity_cons0_4: - Artificial + 600 lambda1 + 600 lambda2 + 375 lambda3
   + 900 lambda4 <= 625
 Capacity_cons0_5: - Artificial + 800 lambda1 + 150 lambda2 + 650 lambda3
   <= 625
 Capacity_cons0_6: - Artificial + 500 lambda4 <= 625
 Capacity_cons1_0: - Artificial + 600 lambda1 + 900 lambda2 + 575 lambda3
   <= 625
 Capacity_cons1_1: - Artificial + 400 lambda1 + 400 lambda2 + 350 lambda3
   + 875 lambda4 <= 625
 Capacity_cons1_2: - Artificial + 400 lambda1 + 425 lambda2 + 400 lambda3
   <= 625
 Capacity_cons1_3: - Artificial + 75 lambda1 + 75 lambda2 + 75 lambda3
   + 325 lambda4 <= 625
 Capacity_cons1_4: - Artificial + 975 lambda1 + 900 lambda4 <= 625
 Capacity_cons1_5: - Artificial + 50 lambda1 + 700 lambda2 + 200 lambda3
   <= 625
 Capacity_cons1_6: - Artificial + 100 lambda1 + 100 lambda2 + 1000 lambda3
   + 500 lambda4 <= 625
 Capacity_cons2_0: - Artificial + 300 lambda1 + 100 lambda3 + 900 lambda4
   <= 625
 Capacity_cons2_1: - Artificial + 750 lambda1 + 100 lambda2 + 650 lambda3
   + 275 lambda4 <= 625
 Capacity_cons2_2: - Artificial + 100 lambda1 + 75 lambda2 + 100 lambda3
   + 500 lambda4 <= 625
 Capacity_cons2_3: - Artificial + 300 lambda1 + 300 lambda2 + 300 lambda3
   + 50 lambda4 <= 625
 Capacity_cons2_4: - Artificial + 225 lambda1 + 1200 lambda2 + 1425 lambda3
   <= 625
 Capacity_cons2_5: - Artificial + 325 lambda1 + 325 lambda2 + 325 lambda3
   + 1175 lambda4 <= 625
 Capacity_cons2_6: - Artificial + 900 lambda1 + 900 lambda2 <= 625
 convexCon: lambda1 + lambda2 + lambda3 + lambda4 = 1
Bounds
End
