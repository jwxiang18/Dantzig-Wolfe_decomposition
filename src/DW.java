import gurobi.*;

import java.util.ArrayList;
import java.util.Arrays;

public class DW {
    Instance instance;
    int iter = 0;
    //凸组合约束的对偶变量，更新子问题目标函数时需要，初始值为1，为了在第一次迭代时，子问题的目标函数能小于1，能有新的检验数为负的列产生
    double dual_convexCons = 1;
    double[][] dual_CapacityCon;
    ArrayList<Double> SP_totalCost = new ArrayList<>();
    GRBModel RMP ,subProblem;
    GRBVar var_Artificial ;
    GRBVar[][][] var_x;
    ArrayList<GRBVar> var_lambda = new ArrayList<>();
    GRBConstr convexCon ;
    GRBConstr[][] CapCons;

    public DW(Instance instance){
        this.instance = instance;
    }

    public void main() throws GRBException {
        initializeModel();

        System.out.println("         *** Dantzig Wolfe Decomposition starts***");
        optimizePhase1();

        updatePhase2();
        optimizePhase2();

        optimizeFinalRMP();
        System.out.println("         *** Dantzig Wolfe Decomposition ends***");
        reportSolution();
    }

    private void initializeModel() throws GRBException {
        GRBEnv env = new GRBEnv();
        RMP = new GRBModel(env);
        subProblem = new GRBModel(env);

        RMP.set(GRB.IntParam.OutputFlag,0);
        subProblem.set(GRB.IntParam.OutputFlag,0);

        // 添加初始的人工变量，使得算法可以进入第一阶段
        var_Artificial = RMP.addVar(0,GRB.INFINITY,0,GRB.CONTINUOUS,"Artificial");

        // 添加初始的人工约束，使得算法可以进入第一阶段 ，即 -1 * var_Artificial <= capacity_ij;
        CapCons = new GRBConstr[instance.NumOrg][instance.NumDes];
        GRBLinExpr expr = new GRBLinExpr();
        expr.addTerm(-1,var_Artificial);
        for (int i = 0; i < instance.NumOrg; i++) {
            for (int j = 0; j < instance.NumDes; j++) {
                CapCons[i][j] = RMP.addConstr(expr , GRB.LESS_EQUAL ,  instance.Capacity[i][j] , "Capacity_cons" + i + "_" + j);
            }
        }

        // 添加凸组合约束
        convexCon = RMP.addConstr(var_Artificial , GRB.EQUAL , 1 , "convexCon");

        // 添加子问题决策变量 x_ijk
        var_x = new GRBVar[instance.NumOrg][instance.NumDes][instance.NumProduct];
        for (int i = 0; i < instance.NumOrg; i++) {
            for (int j = 0; j < instance.NumDes; j++) {
                for (int k = 0; k < instance.NumProduct; k++) {
                    var_x[i][j][k] = subProblem.addVar(0,GRB.INFINITY,0,GRB.CONTINUOUS,"x" + i + "_" + j + "_" + k);
                }
            }
        }

        for (int i = 0; i < instance.NumOrg; i++) {
            for (int j = 0; j < instance.NumProduct; j++) {
                GRBLinExpr expr_artificial = new GRBLinExpr();
                for (int k = 0; k < instance.NumDes; k++) {
                    expr_artificial.addTerm(1,var_x[i][k][j]);
                }
                subProblem.addConstr(expr_artificial , GRB.EQUAL , instance.Supply[i][j] , "Supply_cons" + i + "_" + j);
            }
        }

        for (int i = 0; i < instance.NumDes; i++) {
            for (int j = 0; j < instance.NumProduct; j++) {
                GRBLinExpr expr1 = new GRBLinExpr();
                for (int k = 0; k < instance.NumOrg; k++) {
                    expr1.addTerm(1,var_x[k][i][j]);
                }
                subProblem.addConstr(expr1 , GRB.EQUAL , instance.Demand[i][j] , "Demand_cons" + i + "_" + j);
            }
        }
    }
    private void optimizePhase1() throws GRBException {

        dual_CapacityCon = new double[instance.NumOrg][instance.NumDes];
        for (int i = 0; i < instance.NumOrg; i++) {
            Arrays.fill(dual_CapacityCon[i] , 0);
        }

        // 设置主问题和子问题的目标函数
        GRBLinExpr obj_master_phase_1 = new GRBLinExpr();
        obj_master_phase_1.addTerm(1,var_Artificial);
        RMP.setObjective(obj_master_phase_1 , GRB.MINIMIZE);

        GRBLinExpr obj_sub_phase_1 = updateObj_sub_phase_1();

        subProblem.setObjective(obj_sub_phase_1 , GRB.MINIMIZE);

        // 为保证最初的模型可行，我们将最初的凸组合约束设置为空 NULL
        RMP.chgCoeff(convexCon , var_Artificial , 0);

        System.out.println("         *** Phase 1 starts***");
        while (true){
            System.out.println("         *** Iteration " + iter + " starts***");

            subProblem.optimize();
            subProblem.write("subProblem.lp");
            if(subProblem.get(GRB.DoubleAttr.ObjVal) >= -1e-6){
                System.out.println("无新的列生成，无检验数为负");
                break;
            }else{
                iter++;
                double SP_totalCost_temp = getTotalCost_subProblem();
                SP_totalCost.add(SP_totalCost_temp);

                GRBColumn col = getGrbColumn();

                var_lambda.add(RMP.addVar(0,GRB.INFINITY,0,GRB.CONTINUOUS,col,"lambda" + iter));

                RMP.optimize();
                RMP.write("RMP.lp");
                if(RMP.get(GRB.DoubleAttr.ObjVal) <= 1e-6){
                    System.out.println("------------------Phase 1 ends------------------");
                    break;
                }else{
                    updateDualValue();

                    obj_sub_phase_1 = updateObj_sub_phase_1();
                    subProblem.setObjective(obj_sub_phase_1 , GRB.MINIMIZE);
                }
            }
        }
    }

    private GRBColumn getGrbColumn() throws GRBException {
        GRBColumn col = new GRBColumn();
        for (int i = 0; i < instance.NumOrg; i++) {
            for (int j = 0; j < instance.NumDes; j++) {
                double flow_i_j = 0;
                for (int k = 0; k < instance.NumProduct; k++) {
                    flow_i_j += var_x[i][j][k].get(GRB.DoubleAttr.X);
                }
                col.addTerm(flow_i_j ,CapCons [i][j]);
            }
        }
        // 注意有凸组合约束对应的列
        col.addTerm(1 , convexCon);
        return col;
    }

    private double getTotalCost_subProblem() throws GRBException {
        double totalCost = 0;
        for (int i = 0; i < instance.NumOrg; i++) {
            for (int j = 0; j < instance.NumDes; j++) {
                for (int k = 0; k < instance.NumProduct; k++) {
                    totalCost += var_x[i][j][k].get(GRB.DoubleAttr.X) * instance.Cost[i][j][k];
                }
            }
        }
        return totalCost;
    }

    private void updateDualValue() throws GRBException {
        for (int i = 0; i < instance.NumOrg; i++) {
            for (int j = 0; j < instance.NumDes; j++) {
                dual_CapacityCon[i][j] = CapCons[i][j].get(GRB.DoubleAttr.Pi);
            }
        }
        dual_convexCons = convexCon.get(GRB.DoubleAttr.Pi);
    }

    private GRBLinExpr updateObj_sub_phase_1() {
        GRBLinExpr obj_sub_phase_1 = new GRBLinExpr();
        for (int i = 0; i < instance.NumOrg; i++) {
            for (int j = 0; j < instance.NumDes; j++) {
                for (int k = 0; k < instance.NumProduct; k++)
                    obj_sub_phase_1.addTerm(-1 * dual_CapacityCon[i][j] , var_x[i][j][k]);
            }
        }
        obj_sub_phase_1.addConstant(-1*dual_convexCons);
        return obj_sub_phase_1;
    }

    private GRBLinExpr updateObj_sub_phase_2() {
        GRBLinExpr obj_sub_phase_2 = new GRBLinExpr();
        for (int i = 0; i < instance.NumOrg; i++) {
            for (int j = 0; j < instance.NumDes; j++) {
                for (int k = 0; k < instance.NumProduct; k++) {
                    obj_sub_phase_2.addTerm((instance.Cost[i][j][k] - dual_CapacityCon[i][j]) , var_x[i][j][k]);
                }
            }
        }
        obj_sub_phase_2.addConstant(-1*dual_convexCons);
        return obj_sub_phase_2;
    }

    private void updatePhase2() throws GRBException {
        System.out.println("         *** Phase 2 starts***");

        GRBLinExpr obj_master_phase_2 = new GRBLinExpr();
        for (int i = 0; i < var_lambda.size(); i++) {
            obj_master_phase_2.addTerm(SP_totalCost.get(i) , var_lambda.get(i));
        }

        RMP.setObjective(obj_master_phase_2 , GRB.MINIMIZE);

        var_Artificial.set(GRB.DoubleAttr.LB , 0);
        var_Artificial.set(GRB.DoubleAttr.UB , 0);

        RMP.optimize();

        updateDualValue();

        GRBLinExpr obj_sub_phase_2 = updateObj_sub_phase_2();

        subProblem.setObjective(obj_sub_phase_2 , GRB.MINIMIZE);
        iter++;
    }

    private void optimizePhase2() throws GRBException {
        System.out.println("------start optimizePhase2------");
        while(true){
            System.out.println("         *** Iteration " + iter + " starts***");

            subProblem.optimize();
            if(subProblem.get(GRB.DoubleAttr.ObjVal) >= -1e-6) {
                System.out.println("无新的列生成，无检验数为负");
                break;
            }else{
                iter++;
                double totalCost_subProblem = getTotalCost_subProblem();

                GRBColumn col = getGrbColumn();
                var_lambda.add(RMP.addVar(0,GRB.INFINITY,totalCost_subProblem,GRB.CONTINUOUS,col,"lambda" + iter));
                RMP.optimize();

                updateDualValue();

                GRBLinExpr obj_sub_phase_2 = updateObj_sub_phase_2();
                subProblem.setObjective(obj_sub_phase_2 , GRB.MINIMIZE);
            }
        }
    }

    private void optimizeFinalRMP() throws GRBException {
        // 主问题中没有存储每个方案具体的流量，因此将主问题的最优解带入子问题中，求出每个方案的流量
        double [][] opt_x = new double[instance.NumOrg][instance.NumDes];
        for (int i = 0; i < instance.NumOrg; i++) {
            Arrays.fill(opt_x[i],0);
        }
        System.out.println("------start optimizeFinalRMP------");
        for (int i = 0; i < instance.NumOrg; i++) {
            for (int j = 0; j < instance.NumDes; j++) {
                opt_x[i][j] = instance.Capacity[i][j] + var_Artificial.get(GRB.DoubleAttr.X) - CapCons[i][j].get(GRB.DoubleAttr.Slack);
            }
        }
        GRBLinExpr obj_sub_final = new GRBLinExpr();
        for (int i = 0; i < instance.NumOrg; i++) {
            for (int j = 0; j < instance.NumDes; j++) {
                for (int k = 0; k < instance.NumProduct; k++) {
                    obj_sub_final.addTerm(instance.Cost[i][j][k] , var_x[i][j][k]);
                }
            }
        }
        subProblem.setObjective(obj_sub_final , GRB.MINIMIZE);

        // 最后的子问题添加约束，从i到j的流量之和等于主问题中的流量
        for (int i = 0; i < instance.NumOrg ; i++) {
            for (int j = 0; j < instance.NumDes; j++) {
                GRBLinExpr sumK = new GRBLinExpr();
                for (int k = 0; k < instance.NumProduct; k++) {
                    sumK.addTerm(1 , var_x[i][j][k]);
                }
                subProblem.addConstr(sumK , GRB.EQUAL , opt_x[i][j] , "CapacityCon" + i + j);
            }
        }
        subProblem.optimize();
    }
    private void reportSolution() throws GRBException {
        System.out.println("------start reportSolution------");
        System.out.println("最优解为：" + subProblem.get(GRB.DoubleAttr.ObjVal));
        System.out.println("最优解对应的方案为：");
        for (int i = 0; i < instance.NumOrg; i++) {
            for (int j = 0; j < instance.NumDes; j++) {
                for (int k = 0; k < instance.NumProduct; k++) {
                    if(var_x[i][j][k].get(GRB.DoubleAttr.X) > 1e-6){
                        System.out.println("从" + i + "到" + j + "的" + k + "的流量为：" + Math.round(var_x[i][j][k].get(GRB.DoubleAttr.X)*100)/100.0);
                    }
                }
            }
        }
    }
}
