from Data import Read_data
from gurobipy import *
import copy
class Solve:
    def __init__(self):
        #input
        self.multiplier=1
        data=Read_data(self.multiplier)
        self.node_list, self.link_list, self.candidate_link_list, self.OD_pair_list, \
        self.g_number_of_links, self.g_number_of_nodes, self.g_number_of_candidate_links, self.g_number_of_ODs=data.read_candidate_links()
        self.reliability=1
        self.construction_budget=0
        self.construction_cost=1

        #Parameter: CG time
        self.iteration_times=1000

        #columns
        self.solutions_of_routing_subproblem = [[] for i in range(self.g_number_of_ODs)] #columns
        self.candidate_link_flag_of_routing_subproblem = [[] for i in range(self.g_number_of_ODs)]  # columns
        self.primal_cost_of_routing_subproblem = [[] for i in range(self.g_number_of_ODs)]
        self.solutions_of_KS_subproblem = []  # columns

        #output
        self.reduced_cost_of_routing_subproblem = [[] for i in range(self.g_number_of_ODs)]
        self.reduced_cost_of_KS_subproblem=[]
        self.optimal_routing_result=[]
        self.optimal_construction_result=[]
        self.record_multiplier_miu = []  # dual variables
        self.global_LB = None
        self.global_UB = 1000000
        self.obj_of_RMP=[]

    def g_solving_RNDP_by_CG(self):
        print("Solving...")
        # solution=[0] * self.g_number_of_candidate_links
        solution, obj = self.g_solving_KS()
        # construction
        self.solutions_of_KS_subproblem.append(solution)

        #Initial columns
        self.g_produce_UB_based_on_KS(solution)

        # step 1: lower bound generation (CG)
        for i in range(self.iteration_times):
            print("iteration:{}".format(i+1))
            Flag=0
            self.record_multiplier_miu.append([])
            # solve the relaxed master problem by LP
            obj_of_RMP, pi_for_construction, pi_list_for_flow_balance = self.g_solving_RMP_by_LP(i)
            self.obj_of_RMP.append(obj_of_RMP)
            # solve the routing subproblems and add columns
            for od_index in range(self.g_number_of_ODs):
                od_pair = self.OD_pair_list[od_index]
                node_seq,candidate_travel_flag,obj,primal_cost = self.g_solving_RSP_version_II(od_index, od_pair, 1)
                self.solutions_of_routing_subproblem[od_index].append(node_seq)
                self.candidate_link_flag_of_routing_subproblem[od_index].append(candidate_travel_flag)
                reduced_cost=round(obj+pi_list_for_flow_balance[od_index],3)
                self.reduced_cost_of_routing_subproblem[od_index].append(reduced_cost)
                self.primal_cost_of_routing_subproblem[od_index].append(primal_cost)
                if reduced_cost<0:
                    Flag=1

            # solve the Knapsack problem and add column
            solution, obj = self.g_solving_KS()
            self.solutions_of_KS_subproblem.append(solution)
            reduced_cost=round(obj+pi_for_construction,3)
            self.reduced_cost_of_KS_subproblem.append(reduced_cost)
            if reduced_cost < 0:
                Flag = 1

            # self.g_produce_UB_based_on_KS(solution)
            # check the reduced costs

            if Flag == 0:
                self.global_LB=obj_of_RMP
                break

        #step 2: Generate an upper bound and output the optimal solution
        obj_of_RMP=self.g_solving_RMP_by_IP()
        if self.global_UB>obj_of_RMP:
            self.global_UB=obj_of_RMP
        print("Lb:{}".format(self.global_LB))
        print("Ub:{}".format(self.global_UB))

    def g_solving_RMP_by_LP(self,i):
        self.RMP=Model("Relaxed master problem")
        self.RMP.setParam('OutputFlag', 0)
        # variable
        for k in range(self.g_number_of_ODs):
            columns_for_OD_k=self.solutions_of_routing_subproblem[k]
            for l in range(len(columns_for_OD_k)):
                self.RMP.addVar(lb=0,ub=1,obj=self.primal_cost_of_routing_subproblem[k][l],vtype=GRB.CONTINUOUS,name="x_{}_{}".format(k,l))

        for m in range(len(self.solutions_of_KS_subproblem)):
            self.RMP.addVar(lb=0,ub=1,obj=0,name="y_{}".format(m))

        self.RMP.update()
        # constraint I: we should select a column for each OD pair
        for k in range(self.g_number_of_ODs):
            columns_for_OD_k = self.solutions_of_routing_subproblem[k]
            expr=LinExpr()
            for l in range(len(columns_for_OD_k)):
                name=self.RMP.getVarByName("x_{}_{}".format(k,l))
                expr.addTerms(1,name)
            self.RMP.addConstr(expr,GRB.EQUAL,1,name="flow balance for OD-{}".format(k))

        # constraint II: select a construction
        expr=LinExpr()
        for m in range(len(self.solutions_of_KS_subproblem)):
            name=self.RMP.getVarByName("y_{}".format(m))
            expr.addTerms(1,name)
        self.RMP.addConstr(expr, GRB.EQUAL, 1, name="flow balance for construction")

        #constraints III: coupling constraints
        for k in range(self.g_number_of_ODs):
            candidate_link_flag_for_OD_k = self.candidate_link_flag_of_routing_subproblem[k]
            for link_id in range(self.g_number_of_candidate_links):
                expr=LinExpr()
                #x
                for l in range(len(candidate_link_flag_for_OD_k)):
                    value=candidate_link_flag_for_OD_k[l][link_id]
                    name=self.RMP.getVarByName("x_{}_{}".format(k,l))
                    expr.addTerms(value,name)
                #y
                for m in range(len(self.solutions_of_KS_subproblem)):
                    solution=self.solutions_of_KS_subproblem[m]
                    value=solution[link_id]*-1
                    name = self.RMP.getVarByName("y_{}".format(m))
                    expr.addTerms(value,name)
                self.RMP.addConstr(expr,GRB.LESS_EQUAL,0,name="Coupling_{}_{}".format(link_id,k))

        self.RMP.optimize()
        self.RMP.write("RMP.lp")
        obj_of_RMP=self.RMP.objval

        # update and record the dual variables!
        for link_id in range(self.g_number_of_candidate_links):
            multiplier_link = []
            for k in range(self.g_number_of_ODs):
                constr = self.RMP.getConstrByName("Coupling_{}_{}".format(link_id, k))
                Pi = -1*constr.Pi #>0
                multiplier_link.append(Pi)
            self.record_multiplier_miu[i].append(multiplier_link)
            self.candidate_link_list[link_id].base_profit_for_lagrangian=multiplier_link

        #should we record its solutions? unnecessary
        constr = self.RMP.getConstrByName("flow balance for construction")
        pi_for_construction = -1*constr.Pi
        pi_list_for_flow_balance=[]
        for k in range(self.g_number_of_ODs):
            constr = self.RMP.getConstrByName("flow balance for OD-{}".format(k))
            pi_list_for_flow_balance.append(-1*constr.Pi)

        return obj_of_RMP,pi_for_construction,pi_list_for_flow_balance

    def g_solving_RMP_by_IP(self):
        self.RMP=Model("Relaxed master problem")
        self.RMP.setParam('OutputFlag', 0)
        # variable
        for k in range(self.g_number_of_ODs):
            columns_for_OD_k=self.solutions_of_routing_subproblem[k]
            for l in range(len(columns_for_OD_k)):
                self.RMP.addVar(obj=self.primal_cost_of_routing_subproblem[k][l],vtype=GRB.BINARY,name="x_{}_{}".format(k,l))

        for m in range(len(self.solutions_of_KS_subproblem)):
            self.RMP.addVar(obj=0,vtype=GRB.BINARY,name="y_{}".format(m))

        self.RMP.update()
        # constraint I: we should select a column for each OD pair
        for k in range(self.g_number_of_ODs):
            columns_for_OD_k = self.solutions_of_routing_subproblem[k]
            expr=LinExpr()
            for l in range(len(columns_for_OD_k)):
                name=self.RMP.getVarByName("x_{}_{}".format(k,l))
                expr.addTerms(1,name)
            self.RMP.addConstr(expr,GRB.EQUAL,1,name="flow balance for OD-{}".format(k))

        # constraint II: select a construction
        expr=LinExpr()
        for m in range(len(self.solutions_of_KS_subproblem)):
            name=self.RMP.getVarByName("y_{}".format(m))
            expr.addTerms(1,name)
        self.RMP.addConstr(expr, GRB.EQUAL, 1, name="flow balance for construction")

        #constraints III: coupling constraints
        for k in range(self.g_number_of_ODs):
            candidate_link_flag_for_OD_k = self.candidate_link_flag_of_routing_subproblem[k]
            for link_id in range(self.g_number_of_candidate_links):
                expr=LinExpr()
                #x
                for l in range(len(candidate_link_flag_for_OD_k)):
                    value=candidate_link_flag_for_OD_k[l][link_id]
                    name=self.RMP.getVarByName("x_{}_{}".format(k,l))
                    expr.addTerms(value,name)
                #y
                for m in range(len(self.solutions_of_KS_subproblem)):
                    solution=self.solutions_of_KS_subproblem[m]
                    value=solution[link_id]*-1
                    name = self.RMP.getVarByName("y_{}".format(m))
                    expr.addTerms(value,name)
                self.RMP.addConstr(expr,GRB.LESS_EQUAL,0,name="Coupling_{}_{}".format(link_id,k))

        self.RMP.optimize()
        obj_of_RMP=self.RMP.objval
        solution=self.RMP.getVars()
        #recover the optimal solution
        index=0
        for k in range(self.g_number_of_ODs):
            columns_for_OD_k=self.solutions_of_routing_subproblem[k]
            for l in range(len(columns_for_OD_k)):
                if round(solution[index].x)==1:
                    self.optimal_routing_result.append(copy.copy(columns_for_OD_k[l]))
                index+=1
        index=-1
        for m in range(len(self.solutions_of_KS_subproblem)):
            if round(solution[index].x) == 1:
                self.optimal_construction_result=copy.copy(self.solutions_of_KS_subproblem[index])
            index-=1
        return obj_of_RMP

    def g_produce_UB_based_on_KS(self,solution):
        for link_id in range(self.g_number_of_candidate_links):
            value = solution[link_id]
            link = self.candidate_link_list[link_id]
            if value == 1:
                link.construction_Flag = 1
            else:
                link.construction_Flag = 0

        local_UB=0
        for od_index in range(self.g_number_of_ODs):
            od_pair = self.OD_pair_list[od_index]
            node_seq, candidate_travel_flag, reduced_cost, primal_cost = self.g_solving_RSP_version_II(od_index,od_pair, 2)
            self.solutions_of_routing_subproblem[od_index].append(node_seq)
            self.candidate_link_flag_of_routing_subproblem[od_index].append(candidate_travel_flag)
            # self.reduced_cost_of_routing_subproblem[od_index].append(reduced_cost)
            self.primal_cost_of_routing_subproblem[od_index].append(primal_cost)
            local_UB+=primal_cost

        if self.global_UB>local_UB:
            self.global_UB=local_UB

    def g_solving_RSP_version_II(self,od_index,od_pair,Flag):
        self.RSP=Model("RSP")
        self.RSP.setParam('OutputFlag', 0)
        expr=LinExpr()
        for link in self.link_list:
            name = "x_{}_{}".format(link.from_node_id, link.to_node_id)
            name = self.RSP.addVar(vtype=GRB.BINARY, name=name)

            value=link.travel_time_mean

            if Flag==1 and link.link_type == 1:
                value+=link.base_profit_for_lagrangian[od_index]

            # if Flag==2 and link.construction_Flag==1:
            #     value += link.travel_time_mean

            if Flag==2 and link.construction_Flag==0:
                value += 10000

            expr.addTerms(value,name)
        y=self.RSP.addVar(vtype=GRB.CONTINUOUS, name="y", lb=0)
        expr.addTerms(self.reliability,y)
        self.RSP.setObjective(expr,GRB.MINIMIZE)
        self.RSP.update()

        # Flow balance
        for node in self.node_list:
            expr = LinExpr()

            # Flow out
            for outbound_link in node.outbound_links_list:
                name = self.RSP.getVarByName("x_{}_{}".format(outbound_link.from_node_id, outbound_link.to_node_id))
                expr.addTerms(1, name)

            # Flow in
            for inbound_link in node.inbound_links_list:
                name = self.RSP.getVarByName("x_{}_{}".format(inbound_link.from_node_id, inbound_link.to_node_id))
                expr.addTerms(-1, name)

            if node.node_id == od_pair[0]:
                self.RSP.addConstr(expr, GRB.EQUAL, 1, name="Node_{}".format(od_pair[0]))

            if node.node_id == od_pair[1]:
                self.RSP.addConstr(expr, GRB.EQUAL, -1, name="Node_{}".format(od_pair[1]))

            if node.node_id != od_pair[0] and node.node_id != od_pair[1]:
                self.RSP.addConstr(expr, GRB.EQUAL, 0, name="Node_{}".format(node.node_id))

        # variance limit
        expr = LinExpr()
        for link in self.link_list:
            name = self.RSP.getVarByName("x_{}_{}".format(link.from_node_id, link.to_node_id))
            expr.addTerms(link.travel_time_variance, name)

        self.RSP.addConstr(expr, GRB.LESS_EQUAL, y ** 2, name="variance limit")
        self.RSP.update()
        self.RSP.write("RSP.lp")
        self.RSP.optimize()
        obj=self.RSP.objval
        values=self.RSP.getVars()
        node_seq, candidate_travel_flag,primal_cost = self.values_transition(values[:-1], od_pair)

        return node_seq,candidate_travel_flag,obj,primal_cost

    def g_solving_KS(self):
        self.KS=Model("KS")
        self.KS.setParam('OutputFlag', 0)
        expr=LinExpr()
        for link in self.candidate_link_list:
            name = "y_{}_{}".format(link.from_node_id, link.to_node_id)
            name=self.KS.addVar(vtype=GRB.BINARY, name=name)
            value=-sum(link.base_profit_for_lagrangian)#+self.construction_cost
            expr.addTerms(value, name)
        self.KS.setObjective(expr,GRB.MINIMIZE)
        self.KS.update()

        expr = LinExpr()
        for link in self.candidate_link_list:
            name=self.KS.getVarByName("y_{}_{}".format(link.from_node_id, link.to_node_id))
            expr.addTerms(1,name)
        self.KS.addConstr(expr,GRB.LESS_EQUAL,self.construction_budget)
        self.KS.optimize()
        obj=self.KS.objval
        values=self.KS.getVars()
        solution=[]
        for value in values:
            solution.append(value.x)
        return solution,obj

    def values_transition(self,values,od_pair):
        #input: values;output:node seq,use seq
        #path
        path_links = {}
        candidate_travel_flag=[]
        path_mean=0
        path_var=0
        for link in self.link_list:
            link_index = link.link_id
            from_node = link.from_node_id
            to_node = link.to_node_id
            if round(values[link_index].x) == 1:
                path_links[from_node] = to_node
                path_mean+=link.travel_time_mean
                path_var+=link.travel_time_variance
            link_type = link.link_type

            if link_type==1 and round(values[link_index].x) == 1:
                candidate_travel_flag.append(1)
            if link_type==1 and round(values[link_index].x) == 0:
                candidate_travel_flag.append(0)

        node_seq = [od_pair[0]+1]
        current_node = od_pair[0]
        while current_node != od_pair[1]:
            current_node = path_links[current_node]
            node_seq.append(current_node+1)

        primal_cost=path_mean+self.reliability*(path_var)**0.5

        return node_seq,candidate_travel_flag,primal_cost

    def output_results(self,spend_time):
        #LB,UB,GAP,TIME,OTERATION
        with open("optimal_solution.csv","w") as fl:
            fl.write("UB:{}\n".format(self.global_UB))
            fl.write("LB:{}\n".format(self.global_LB))
            gap=(self.global_UB-self.global_LB)/self.global_UB
            fl.write("gap:{}\n".format(gap))
            fl.write("Time:{} seconds\n".format(spend_time))
            fl.write("optimal construction:{}\n".format(self.optimal_construction_result))
            for od_pair in range(self.g_number_of_ODs):
                fl.write("OD_{}:{}\n".format(od_pair,self.optimal_routing_result[od_pair]))

        with open("RMP_obj.csv","w") as fl:
            fl.write("iteration,obj_of_RMP\n")
            for i in range(len(self.obj_of_RMP)):
                fl.write(str(i+1)+","+str(self.obj_of_RMP[i])+"\n")

        with open("dual_values.csv","w") as fl:
            fl.write("i,")
            for k in range(self.g_number_of_ODs):
                for linkid in range(self.g_number_of_candidate_links):
                    fl.write("{}_{},".format(k,linkid))
            fl.write("\n")
            for i in range(len(self.record_multiplier_miu)):
                fl.write(str(i+1)+",")
                for k in range(self.g_number_of_ODs):
                    for linkid in range(self.g_number_of_candidate_links):
                        multiplier=round(self.record_multiplier_miu[i][linkid][k],3)
                        fl.write(str(multiplier)+",")
                fl.write("\n")

        with open("Reduced_cost.csv","w") as fl:
            fl.write("ite,KS,")
            for od_pair in range(self.g_number_of_ODs):
                fl.write("RSP_{},".format(od_pair))
            fl.write("\n")

            for i in range(len(self.reduced_cost_of_KS_subproblem)):
                fl.write(str(i+1)+",")
                fl.write(str(self.reduced_cost_of_KS_subproblem[i])+",")
                for od_pair in range(self.g_number_of_ODs):
                    fl.write(str(self.reduced_cost_of_routing_subproblem[od_pair][i])+",")
                fl.write("\n")