package Algo;

import base.Customer;
import base.Route;
import lombok.Data;
import param.Conf;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * @description 染色体类，可以通过split函数和Solution类进行转化
 * 转化到 solution后 solution 中的方法得到 irank ,idis ,也可以在此返回到 适应度
 */
@Data
public class Chromosome {
    final int INF = 999999;
    public int index = -1;//作为唯一标识的话，中间变量的ID怎么办？
    public Conf userConf;

    public ArrayList<Integer> cur_list;

    public double fitness = 0;
    public double chrProfit;
    public double chrDistance;
    public double chrUsedVehicle;
    public int irank = 0;
    public double idis = 0;
    //public double crowd1 = 0;
    //public double crowd2 = 0;

    public int np;
    public ArrayList<Chromosome> sp_list;

    public int preIndex;
    public Customer depot;

    public ArrayList<double[]> Mts = new ArrayList<>();
    public ArrayList<Double> Prs = new ArrayList<>();

    FileWriter fw = new FileWriter("./dataset/chromosomeProbabilityLog.txt");
    BufferedWriter bw = new BufferedWriter(fw);
    PrintWriter pw = new PrintWriter(bw);

    //构造函数，随机生成一个初始解
    public Chromosome(Conf userConf) throws IOException {
        this.userConf = userConf;
        cur_list = new ArrayList<>();
        for (int i = 1; i <= userConf.N; i++) {
            this.cur_list.add(i);
        }
        Collections.shuffle(this.cur_list);
        //todo 随机删除几个点变成N+1,N+2...;
        unvisited(this.cur_list,new Random().nextInt(5));
        cur_list.add(0, 0);
        //System.out.println(cur_list.toString());
        depot = userConf.customers[0];
        //初始解多shuffle 几次再改进
    }

    private void unvisited(ArrayList<Integer> cur_list, int cnt) {
        Random r = new Random();
        //int cnt = r.nextInt(userConf.N / 2);
        //cnt = 2;
        for (int i = 0; i < cnt; i++) {
            int j = r.nextInt(cur_list.size());
            cur_list.set(j, userConf.N + 1 + i);
        }
    }

    //使用分割函数：跑一遍bellman-ford算法获得最优分割
    //实际上转化为从开始点到结束点的最短路划分问题
    //将编码转换为solution
    public Solution toSolution() {
        Solution solution = new Solution(userConf);
        double[] V = new double[userConf.N + 1 + userConf.newN];//距离数组
        int[] P = new int[userConf.N + 1 + userConf.newN]; // 储存了连接该点的上一点
        int j; // 循环的标识
        int cost;// 当前的花费
        int dist; //
        double time;// 当前的时间
        for (int i = 1; i <= userConf.N + userConf.newN; i++)//到达的点的最少花费
            V[i] = INF;
        //int index = 0;
        for (int i = 1; i <= userConf.N + userConf.newN; i++) {
            //P[i] = this.cur_list.get(i);//最开始所有点都没连上
            P[i] = i;
            //System.out.println("P: " + i + "," + P[i]);
        }
        for (int i = 1; i <= userConf.N; i++) {
            cost = 0;//demand
            time = 0;//travel time cost
            j = i;
            // 重新放一遍
            Mts = new ArrayList<>();
            Prs = new ArrayList<>();
            preIndex = 0;
            boolean flag = true;
            while (true) {
                //行程的开始
                if (i == j) {
                    time += Math.max(userConf.customers[cur_list.get(j)].r_time, userConf.dis_matriax[0][cur_list.get(j)]);
                    time += userConf.customers[cur_list.get(j)].s_time; // 服务时间
                    time += userConf.dis_matriax[cur_list.get(j)][0];
                    cost += userConf.customers[cur_list.get(j)].demand;
                    // add 0
                    depot = userConf.customers[0];
                    double[] Mt0_ = new double[depot.d_time * 2 + 1];
                    for (int i1 = 0; i1 < Mt0_.length; i1++) {
                        Mt0_[i1] = 1.0;
                    }
                    Mts.add(Mt0_);
                    Prs.add(Mt0_[depot.d_time]);
                    int preCustomerIndex = 0;
                    int curCustomerIndex = cur_list.get(j);
                    flag = this.check_time(Mts, preCustomerIndex, curCustomerIndex);
                    pw.print("\n Start pre: " + this.cur_list.get(j - 1) + " cur: " + cur_list.get(j)
                            + ",Pr= " + Prs.toString() + " flag= " + flag + " preINdex" + preIndex);
                    pw.flush();
                } else {
                    double next_time = time - userConf.dis_matriax[cur_list.get(j - 1)][0] + userConf.dis_matriax[cur_list.get(j)][cur_list.get(j - 1)];//到达下一个的时间
                    int preCustomerIndex = cur_list.get(j - 1);
                    int curCustomerIndex = cur_list.get(j);
                    flag = this.check_time(Mts, preCustomerIndex, curCustomerIndex);
                    pw.print("\n ELSE pre: " + this.cur_list.get(j - 1) + " cur: " + cur_list.get(j)
                            + ",Pr= " + Prs.toString() + " flag= " + flag + " preINdex" + preIndex);
                    pw.flush();
                    //if (next_time > userConf.customers[cur_list.get(j)].d_time || !flag)
                    if (!flag)
                        break;//
                    time = Math.max(next_time, userConf.customers[cur_list.get(j)].r_time);
                    time += userConf.dis_matriax[cur_list.get(j)][0];
                    cost += userConf.customers[cur_list.get(j)].demand;
                }
                //假如满足容量约束和时间约束
                if (cost <= userConf.Cap && time <= userConf.maxrouteLength && flag) { //---容量约束未被违反的情况下 可以改成每段route-限制长度
                    if (V[cur_list.get(j)] > V[cur_list.get(i - 1)] + time) {
                        V[cur_list.get(j)] = V[cur_list.get(i - 1)] + time;//不断更新当前最短路
                        P[cur_list.get(j)] = this.cur_list.get(i - 1);
                        pw.print("\n pre: " + this.cur_list.get(i - 1) + " cur: " + cur_list.get(j)
                                + ",Pr= " + Prs.toString() + " flag= " + flag);
                        pw.flush();
                    }
                    j++;
                }
                if (j > userConf.N || time >= userConf.customers[0].d_time || cost >= userConf.Cap)
                    break;
            }
        }
        //System.out.println("After: ");
        Route route = new Route();
        for (int i = 1; i <= userConf.N; i++) {
            //System.out.println("P: " + i + "," + P[i]);
        }
        int tmp = P[cur_list.get(userConf.N)];
        int i = userConf.N;
        // 将分割过的重新组成Solution
        while (i > 0) {
            if (P[cur_list.get(i)] == tmp)
                route.cus_list.add(cur_list.get(i));
            else {
                tmp = P[cur_list.get(i)];
                route.getVehicle();
                route.getValue(userConf);
                route.getProfit(userConf);
                Collections.reverse(route.cus_list);
                solution.rou_list.add(route);
                route = new Route();
                route.cus_list.add(cur_list.get(i));
            }
            i--;
        }
        if (route.cus_list.size() != 0) {
            Collections.reverse(route.cus_list);
            route.getVehicle();
            route.getValue(userConf);
            route.getProfit(userConf);
            solution.rou_list.add(route);
        }
        return solution;//check
    }

    public boolean check_time(ArrayList<double[]> Mts,
                              int preCustomerIndex, int curCustomerIndex) {

        Customer curC = userConf.customers[curCustomerIndex];
        double[] Mtj_ = new double[depot.d_time * 2 + 1];
        for (int z = 0; z < curC.r_time; z++) {
            Mtj_[z] = 0;
        }
        double[] Mti_ = Mts.get(preIndex);//todo
        for (int z = curC.r_time; z <= curC.d_time; z++) {
            for (Integer k : userConf.mik[preCustomerIndex]) {
                int idx1 = (int) (z - k - userConf.dis_matriax[curCustomerIndex][preCustomerIndex]);
                if (idx1 < 0) {
                    Mtj_[z] += 0;
                } else {
                    //System.out.println("curC  " + curCustomerIndex + ",k= " + k + ",idx1= " + idx1);
                    Mtj_[z] += userConf.ms[preCustomerIndex][k] * Mti_[idx1];
                }
            }
        }
        for (int z = curC.d_time + 1; z <= 2 * depot.d_time; z++) {
            Mtj_[z] = Mtj_[curC.d_time];
        }
        double Pj = Mtj_[curC.d_time];
        Mts.add(Mtj_);
        Prs.add(Pj);
        preIndex += 1;
        return (Pj >= userConf.alpha);

    }

    public void toSolutionS() {
    }

    public Chromosome copy() throws IOException {
        Chromosome chromosome = new Chromosome(this.userConf);
        chromosome.cur_list.clear();
        chromosome.cur_list.addAll(this.cur_list);
        chromosome.index = this.index;
        chromosome.fitness = this.fitness;
        return chromosome;
    }

    // 重新返回--->设置fitness
    //public void setFitness() {
    //this.fitness = this.toSolution().getFitness();
    //}

    public void setFitness1() {
        this.fitness = this.toSolution().getFitness();
        this.chrDistance = this.toSolution().solutionDistance;
        this.chrProfit = this.toSolution().solutionProfit;
        this.chrUsedVehicle = this.toSolution().solutionUsedVehicle;
        //System.out.println("idx,obj= " + this.index + "," + this.fitness);
    }

    public void print() {
        this.toSolution().print();
        System.out.println("idis,irank = " + this.irank);
    }
}
