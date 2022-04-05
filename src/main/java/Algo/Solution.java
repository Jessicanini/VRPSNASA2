package Algo;

import base.Route;
import param.Conf;
import util.utils;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

/**
 * @description Solution类，构成了问题的一个解,提供了转到chromosome的函数
 * chromosome 转到solution 中间类 为 route
 */
public class Solution {
    public Conf UserConf;

    public ArrayList<Route> rou_list = new ArrayList<>();

    public double fitness; // 适应度
    public double solutionProfit;
    public double solutionDistance;
    public double solutionUsedVehicle;

    public Solution(Conf UserConf) {
        this.UserConf = UserConf;
    }

    // fitness ---> route的dis, profit 车辆数目
    public double getFitness() {
        this.fitness = 0;
        this.solutionProfit = 0;
        this.solutionDistance = 0;
        this.solutionUsedVehicle = 0;

        for (Route route : rou_list) {
            route.getVehicle();
            route.getValue(UserConf);
            route.getProfit(UserConf);
            this.fitness += UserConf.a1 * route.value - UserConf.a2 * route.profit;
            this.solutionProfit += route.profit;
            this.solutionDistance += route.value;
            this.solutionUsedVehicle += route.usedVehicle;
        }
        return this.fitness;
    }

    //todo
    public boolean isParetoDominate(Solution other) {
        this.getFitness();
        other.getFitness();
        double thisP = this.solutionProfit;
        double thisD = this.solutionDistance;
        double otherP = other.solutionProfit;
        double otherD = other.solutionDistance;
        //全部优于，至少1个没有等于
        boolean isdominated = (utils.greatereq(thisP, otherP) && utils.greatereq(otherD, thisD)) &&
                (utils.greater(thisP, otherP) || utils.greater(otherD, thisD));
        //System.out.println("thisP,otherP,thisD,otherD = " + thisP + "," + otherP + "," + thisD + "," + otherD);
        return isdominated;
    }

    //把solution转化为chromosome
    public Chromosome tochromosome() throws IOException {
        Chromosome chromosome = new Chromosome(UserConf);
        for (Route route : rou_list) {
            chromosome.cur_list.addAll(route.cus_list);
        }
        return chromosome;
    }

    //输出解 --- 处理没有访问的顾客的情况
    public void print() {
        int i = 1;
        ArrayList<Integer> unvisited = new ArrayList<>();
        int[] tmp = new int[UserConf.N + 1 + UserConf.newN];
        for (int j = 1; j < UserConf.N + 1 + UserConf.newN; j++) {
            tmp[j] = 0;
        }
        for (Route route : this.rou_list) {
            System.out.print("Route " + i + ": 0-");
            for (int j : route.cus_list) {
                System.out.print(j + "-");
                tmp[j]++;
            }
            System.out.print("0");
            System.out.print(" 可行性检验: " + route.check(this.UserConf));
            System.out.print(" 路径长度: " + route.value + " 单条路径利润: " + route.profit);
            System.out.println();
            i++;
        }
        for (int k = 1; k <= UserConf.N; k++) {
            if (tmp[k] == 0) {
                unvisited.add(k);
            }
        }
        System.out.println("该解的单目标函数值为" + this.getFitness() + " 收集到的利润 = " +
                this.solutionProfit + " 路线距离 = " + this.solutionDistance
                + " 未访问的顾客:  " + unvisited.toString());
    }

    //输出解到文件 --- 处理没有访问的顾客的情况
    public void write(PrintWriter pw) {
        int i = 1;
        ArrayList<Integer> unvisited = new ArrayList<>();
        int[] tmp = new int[UserConf.N + 1 + UserConf.newN];
        for (int j = 1; j < UserConf.N + 1 + UserConf.newN; j++) {
            tmp[j] = 0;
        }
        for (Route route : this.rou_list) {
            pw.print("[0,");
            for (int j : route.cus_list) {
                pw.print(j + ",");
                tmp[j]++;
            }
            pw.print("0]");
            route.check(this.UserConf);
            //pw.print(" 可行性检验: " + route.check(this.UserConf));
            //pw.print(" 路径长度: " + route.value + " 单条路径利润: " + route.profit);
            pw.print(",");
            i++;
        }
        for (int k = 1; k <= UserConf.N; k++) {
            if (tmp[k] == 0) {
                unvisited.add(k);
            }
        }
        pw.println("该解的单目标函数值为" + this.getFitness() + " 收集到的利润 = " +
                this.solutionProfit + " 路线距离 = " + this.solutionDistance
                + " 未访问的顾客:  " + unvisited.toString());
        pw.flush();
    }

    public void write1(PrintWriter pw, int iter) {
        this.getFitness();
        //pw.println("迭代次数" + iter + " 利润 = " + this.solutionProfit + " 路线距离 = " + this.solutionDistance);
        pw.println(iter + "," + this.solutionProfit + "," + this.solutionDistance);
        pw.flush();
    }


}
