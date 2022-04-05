package base;

import lombok.Data;

/**
 * @description Customer类，储存了顾客的信息
 */
@Data
public class Customer {
    public int x;
    public int y;
    public int demand;
    public int r_time;//开始时间
    public int d_time;//结束时间
    public int s_time;//服务时间
    public double profit = 0.0;
    public int cIndex = -1;

    public Customer() {
    }


}
