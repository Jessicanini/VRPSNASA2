package base;

import param.Conf;

import java.sql.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Collectors;

/**
 * @description Route类，每一个route储存了一辆车所访问的顾客
 */
public class Route {
    public ArrayList<Integer> cus_list = new ArrayList<>();// route列表
    public boolean isfeasible = false;

    public double value;
    public double profit;
    public int usedVehicle;
    public ArrayList<Double> disCost = new ArrayList<>();

    public boolean check(Conf userConf) {
        isfeasible = check_c() && check_stochastictime(userConf);
        //check_t(userConf);
        return isfeasible;
    }

    public boolean check_c() {
        return true;
    }

    //静态时间窗口检查
    public boolean check_t(Conf userConf) {
        double time = 0;
        time += userConf.dis_matriax[0][this.cus_list.get(0)];
        if (time > userConf.customers[cus_list.get(0)].d_time) return false;
        for (int i = 1; i <= this.cus_list.size() - 1; i++) {
            time = Math.max(userConf.customers[this.cus_list.get(i - 1)].r_time, time + userConf.dis_matriax[this.cus_list.get(i - 1)][this.cus_list.get(i)]);
            if (time > userConf.customers[cus_list.get(i)].d_time) return false;
        }
        return true;
    }

    //todo 检查概率条件下的TW
    public boolean check_stochastictime(Conf userConf) {

        Customer depot = userConf.customers[0];

        ArrayList<double[]> Mts = new ArrayList<>();
        ArrayList<Double> Prs = new ArrayList<>();
        double[] Mt0_ = new double[depot.d_time * 2 + 1];
        for (int i = 0; i < Mt0_.length; i++) {
            Mt0_[i] = 1.0;
        }
        Mts.add(Mt0_);
        double p0 = Mt0_[depot.d_time];
        Prs.add(p0);

        int preIndex = 0;
        int curIndex = 0;
        int curCustomerIndex = 0;
        int preCustomerIndex = 0;

        Customer curC = new Customer();
        Customer preC = new Customer();
        ArrayList<Integer> newcus_list = new ArrayList<>();// route列表
        newcus_list = (ArrayList<Integer>) this.cus_list.clone();
        newcus_list.add(0, 0);

        for (int i = 1; i < newcus_list.size(); i++) {
            curIndex = i;
            curCustomerIndex = newcus_list.get(curIndex);
            preCustomerIndex = newcus_list.get(preIndex);
            curC = userConf.customers[curCustomerIndex];
            preC = userConf.customers[preCustomerIndex];
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
            //System.out.println("preC=" + preCustomerIndex + "curC=" + curCustomerIndex + ",Mts= "+Pj);
            //for (int z = curC.r_time; z <= curC.d_time; z++) {
            //System.out.print(z+","+Mtj_[z] + ",");
            //}
            preIndex = curIndex;
        }
        //check Prs
        //System.out.print(Prs.toString());
        int out_time = (int) Prs.stream().filter(c -> c < userConf.alpha).count();
        return (out_time <= 0);
    }

    //获得route的dis
    public double getValue(Conf userConf) {
        this.value = 0;
        value += userConf.dis_matriax[0][cus_list.get(0)];//开始
        value += userConf.dis_matriax[0][cus_list.get(cus_list.size() - 1)];
        if (cus_list.size() > 1) {
            for (int i = 1; i < cus_list.size(); i++) {
                value += userConf.dis_matriax[cus_list.get(i)][cus_list.get(i - 1)];
            }
        }
        return value;
    }

    //获得route的profit
    public double getProfit(Conf userConf) {
        //System.out.println(Arrays.stream(userConf.customers).filter(c -> cus_list.contains(c)).collect(Collectors.toList()).toString());
        profit = Arrays.stream(userConf.customers).filter(c -> cus_list.contains(c.cIndex)).
                mapToDouble(e -> e.getProfit()).sum();
        return profit;
    }

    //获得route的车辆数目（主要保障可行性）
    public double getVehicle() {
        this.usedVehicle = isfeasible ? 1 : 0;
        return this.usedVehicle;
    }

    public void checkDistance(Conf userConf) {
        int pre = 0;
        for (int i = 1; i < this.cus_list.size(); i++) {
            double dis = userConf.dis_matriax[this.cus_list.get(pre)][this.cus_list.get(i)];
            this.disCost.add(dis);
            pre = i;
        }
        double totalDis = this.disCost.stream()
                .collect(Collectors.summingDouble(Double::doubleValue));
        System.out.println("route_dis = " + totalDis + ",dis = " + this.disCost.toString());
    }
}
