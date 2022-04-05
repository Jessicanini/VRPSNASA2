package param;

public class GAParam {
    public int pop_number;// 种群的数量
    public int maxGeneration;//
    public double selectionRate;//
    public double crossoverRate;//
    public double mutateRate;//
    public int ObjM;
    public int maxSize;
    public int maxLocalCnt;
    public int maxInitialCnt;
    public int maxShuffleCnt;

    public GAParam() {

        maxSize = 5;//最多不访问点的个数
        ObjM = 2;//目标函数的个数
        crossoverRate = 0.6;
        mutateRate = 0.001;

        //======
        maxGeneration = 100;
        pop_number = 50;// 种群的数量
        maxLocalCnt = 50;

        maxInitialCnt = 10;
        maxShuffleCnt = 300;
        selectionRate = 0.3;

    }
}
