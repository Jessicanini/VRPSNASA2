package param;

import base.Customer;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

import static java.lang.Math.min;

/**
 * @author minasora
 * @date 2019/10/7 16:16
 * @description 储存问题的公共设置，顾客数，容量约束等,输入输出函数
 */
public class Conf {
    public final int K = 99999;
    public int N; // 顾客数 顾客数量不变
    public int newN;// 虚拟顾客点数
    public Customer[] customers;// 顾客 尝试增加newN个虚拟顾客,表示有点没有被访问

    public String instance_name;// 案例名称
    public double[][] dis_matriax;// 距离矩阵
    //public int[] profit_matriax;// 距离矩阵

    public int Cap; // 容量约束
    public int mina = 9999;
    public int maxb = 0;
    public double maxrouteLength = 99999;// 电动车辆路径长度约束

    public double a1 = 1.0;//距离和利润的系数
    public double a2 = 1.0;

    public double alpha = 0.95;
    public double[][] ms;//=new double[][];// 点 i 服务时长为k  的概率值
    public ArrayList<Integer>[] mik;// 点i 包含的服务时长的数值 1,2,3,...

    public void readInstance(String filename, int a, int c, int b, PrintWriter pw1) throws IOException {
        File file_to_read = new File(filename);
        Scanner cin = new Scanner(file_to_read);
        instance_name = cin.nextLine();
        N = cin.nextInt();
        newN = N;
        Cap = cin.nextInt();
        customers = new Customer[N + 1 + newN];//新建数组
        dis_matriax = new double[N + 1 + newN][N + 1 + newN];
        while (cin.hasNext()) {
            int i = cin.nextInt();
            customers[i] = new Customer();
            customers[i].x = cin.nextInt();
            customers[i].y = cin.nextInt();
            customers[i].demand = cin.nextInt();
            customers[i].r_time = cin.nextInt();
            customers[i].d_time = cin.nextInt();
            customers[i].s_time = cin.nextInt();
            customers[i].profit = cin.nextInt();//todo
            customers[i].cIndex = i;
            mina = Math.min(mina, customers[i].r_time);
            maxb = Math.max(maxb, customers[i].d_time);
        }
        //arificial
        for (int i = 0; i < newN; i++) {
            customers[N + i + 1] = new Customer();
            customers[N + i + 1].x = 0;
            customers[N + i + 1].y = 0;
            customers[N + i + 1].demand = 0;
            customers[N + i + 1].r_time = 0;
            customers[N + i + 1].d_time = maxb * 2;
            customers[N + i + 1].s_time = 0;
            customers[N + i + 1].profit = 0.0;
            customers[N + i + 1].cIndex = N + i + 1;
        }

        // 初始化距离矩阵
        for (int i = 0; i <= N; i++) {
            for (int j = i; j <= N; j++) {
                if (i == j) {
                    dis_matriax[i][j] = 0;
                } else {
                    dis_matriax[i][j] = this.disE(customers[i], customers[j]);
                    dis_matriax[j][i] = dis_matriax[i][j];
                }
            }
        }
        for (int i = 0; i <= N + newN; i++) {
            for (int j = 0; j < newN; j++) {
                dis_matriax[i][N + 1 + j] = 0;
                dis_matriax[N + 1 + j][i] = 0;
            }
        }
        setServiceTime(a, c, b, pw1);
    }

    public void setServiceTime(int a, int c, int b, PrintWriter pw1) {
        // a = 70,c = 90,b=110
        mik = new ArrayList[N + 1 + newN];
        ms = new double[N + 1 + newN][b + 1];
        for (int i = 1; i < N + 1; i++) {
            ArrayList<Integer> tmp = new ArrayList<Integer>();
            for (int j = a; j <= b; j++) {
                tmp.add(j);
            }
            mik[i] = tmp;
            for (int j = a; j <= c; j++) {
                ms[i][j] = 2.0 * (j - a) / ((b - a) * (c - a));
            }
            for (int j = c + 1; j <= b; j++) {
                ms[i][j] = 2.0 * (b - j) / ((b - a) * (b - c));
            }
            pw1.write("\n Node " + i + " , mik= " + mik[i].toString() + " prob ");
            for (int j = a; j <= b; j++) {
                pw1.write(ms[i][j] + ",");
                pw1.flush();
            }
            pw1.write("\n ");
            pw1.flush();
        }
        // 0 -> depot; N+1开始 -> 虚拟顾客
        ArrayList<Integer> tmp = new ArrayList<Integer>();
        tmp.add(0);
        mik[0] = tmp;
        ms[0][0] = 1.0;

        for (int i = N + 1; i < N + 1 + newN; i++) {
            ArrayList<Integer> tmp1 = new ArrayList<Integer>();
            tmp1.add(0);
            mik[i] = tmp1;
            ms[i][0] = 1.0;
        }
    }

    private int disE(Customer a, Customer b) {
        return (int) (Math.sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y)));//返回两个顾客的欧式距离
    }
}
