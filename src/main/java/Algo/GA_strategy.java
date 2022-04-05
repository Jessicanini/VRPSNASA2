package Algo;

import param.Conf;
import param.GAParam;
import util.StdRandom;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;

import static java.util.Collections.max;

/**
 * @author minasora
 * @date 2019/10/7 16:16
 * @description 提供了crossover，mutation，selection等遗传算法使用的函数
 */
public class GA_strategy {
    public int chrIndex = 0;
    public GAParam UserGAParam;
    public Conf UserConf;
    public int pop_number;            //种群大小

    public Chromosome[] parentPop;//父类
    public Chromosome[] childPop;//中间子代 用来取代父类或者和父类合并
    public Chromosome[] tmpPop;
    public Chromosome[] bestPop;//保存每一代里最好的个体
    public List<Chromosome> finalPop;//保存最后一代里非支配的个体

    public double[] selectedProb;
    public double averageFitness;    //平均适应度
    public double[] relativeFitness;    //每一个 个体的适应度
    public double a1_0 = 0.5;
    public double a1_1 = 0.5;
    public int whichlevel = -1;
    public boolean delta = false;

    public Chromosome bestIndividual;//
    public Chromosome currentBest;//

    FileWriter fw, fw1, fw2, fw3;
    BufferedWriter bw, bw1, bw2, bw3;
    PrintWriter pw, pw1, pw2, pw3;


    public GA_strategy(Conf UserConf, GAParam UserGAParam, String prefix) throws IOException {
        this.UserGAParam = UserGAParam;
        this.UserConf = UserConf;
        this.pop_number = UserGAParam.pop_number;
        this.relativeFitness = new double[this.pop_number];

        //fw = new FileWriter("./dataset/" + prefix + "GALog.txt");
        //bw = new BufferedWriter(fw);
        //pw = new PrintWriter(bw);

        //fw1 = new FileWriter("./dataset/" + prefix + "ChromosomeLog.txt");
        //bw1 = new BufferedWriter(fw1);
        //pw1 = new PrintWriter(bw1);

        fw2 = new FileWriter("./dataset/" + prefix + "IterationLog.txt");
        bw2 = new BufferedWriter(fw2);
        pw2 = new PrintWriter(bw2);

        fw3 = new FileWriter("./dataset/" + prefix + "ParetoLog.txt");
        bw3 = new BufferedWriter(fw3);
        pw3 = new PrintWriter(bw3);

        pw3.write("迭代次数 利润 成本\n");
        pw3.flush();
    }

    public Chromosome[] initialize() throws IOException {
        Chromosome[] chromosomes = new Chromosome[pop_number];
        for (int i = 0; i < pop_number; i++) {
            Chromosome tmp = new Chromosome(UserConf);
            Chromosome tmp1 = this.upDateInitialization(tmp);

            chromosomes[i] = tmp1;//包含了initial
            //
            chromosomes[i].setFitness1();
            chromosomes[i].setIndex(chrIndex++);
            System.out.println("Chromosomes Index = " + chrIndex + " Finish ");
            //pw.println("Initial Index = " + chromosomes[i].index + " irank= "
                    //+ chromosomes[i].irank + ",idis= " + chromosomes[i].idis + "," +
                    //chromosomes[i].cur_list.toString());
            //chromosomes[i].toSolution().write(pw);
           // pw.flush();
        }
        return chromosomes;
    }

    public Chromosome upDateInitialization(Chromosome chromosome) throws IOException {
        System.out.println("Chrom IDx= " + chromosome.index);
        double time1 = System.currentTimeMillis();
        chromosome.setFitness1();
        Chromosome oldChromosome = chromosome.copy();
        int iteration = 0;
        System.out.println("打乱前染色体适应度 = " + chromosome.fitness);
        Collections.shuffle(oldChromosome.cur_list.subList(1, this.UserConf.N));
        while (iteration < this.UserGAParam.maxShuffleCnt) {
            oldChromosome.setFitness1();
            if (oldChromosome.fitness < chromosome.fitness) {
                chromosome = oldChromosome.copy();
                //break;
            } else {

            }
            iteration += 1;
            Collections.shuffle(oldChromosome.cur_list.subList(1, this.UserConf.N));
        }
        double time2 = System.currentTimeMillis();
        System.out.println("打乱后染色体适应度 = " + chromosome.fitness +
                " Shuffle iteraion = " + iteration + " cost= " + (time2 - time1) + " ms");

        System.out.println("搜索前染色体适应度 = " + chromosome.fitness);
        Chromosome newchromosome1 = this.neighborhood_move(chromosome);
        int iteration1 = 0;
        while (iteration1 < this.UserGAParam.maxInitialCnt) {
            newchromosome1.setFitness1();
            chromosome.setFitness1();
            //System.out.println("原染色体适应度 = " + chromosome.fitness + ",新 = " + newchromosome1.fitness);
            if (newchromosome1.fitness < chromosome.fitness) {
                chromosome = newchromosome1.copy();
                //break;
            } else {

            }
            iteration1 += 1;
            newchromosome1 = this.neighborhood_move(chromosome);
        }
        double time3 = System.currentTimeMillis();
        System.out.println("搜索后染色体适应度 = " + chromosome.fitness + " Search iteration = " + iteration1 + " cost " + (time3 - time2) + " ms");
        System.out.println();
        return chromosome;
    }

    //选择 --- 二进制锦标赛选择//轮盘赌选择
    public Chromosome[] selection(Chromosome[] chromosomes) {
        Chromosome[] childrens = new Chromosome[pop_number];
        double r = rand();
        if (r < this.UserGAParam.selectionRate) {
            childrens = this.eliteselection(chromosomes);
        } else {
            childrens = this.rouletteWheelSelection(chromosomes);
        }
        //pw.println("selection: " + Arrays.stream(childrens).
               // map(c -> c.index).collect(Collectors.toList()).toString());
        //pw.flush();
        return childrens;
    }

    //OX交叉方式，返回一个新的Chromosome，交换亲代即可获得两个子代
    public Chromosome crossover(Chromosome p_one, Chromosome p_two) throws IOException {
        Chromosome children = new Chromosome(UserConf);
        children.cur_list.clear();//clear 前要保证有
        children.cur_list.add(0);
        Random r = new Random();
        int i = r.nextInt(UserConf.N) + 1;
        int j = r.nextInt(UserConf.N) + 1;
        if (j < i) {
            int mid = i;
            i = j;
            j = mid;
        }
        // 获得 0<= i<= j <= n的随机数
        //System.out.println("fir : "+p_one.cur_list.toString());
        //System.out.println("sec : "+p_two.cur_list.toString());
        //System.out.println("x,y= " + i + "," + j);
        for (int tmp = i; tmp <= j; tmp++)
            children.cur_list.add(p_one.cur_list.get(tmp));//把i-j中的点先转移到mid中
        //System.out.println("child= " +  children.cur_list.toString());
        ArrayList<Integer> parent = new ArrayList<>();
        for (int tmp = j; tmp <= UserConf.N; tmp++)
            parent.add(p_two.cur_list.get(tmp));
        for (int tmp = 1; tmp <= j; tmp++)
            parent.add(p_two.cur_list.get(tmp));
        int tmp = j;
        Boolean if_fir = false;//标记是否到头
        //System.out.println("parent= " +  parent.toString());
        for (int t : parent) {
            if (tmp == UserConf.N) {//假如到了末尾
                tmp = 1;//从头开始
                if_fir = true;
            }
            if (if_fir) {
                if (!children.cur_list.contains(t)) {
                    children.cur_list.add(1, t);
                    tmp++;
                }
            } else {
                if (!children.cur_list.contains(t)) {//不包含的话就加入
                    children.cur_list.add(t);
                    tmp++;
                }
            }
            //System.out.println("child= " +  children.cur_list.toString());
        }
        children.setFitness1();
        children.setIndex(chrIndex++);
        //pw.println("CrossOX: child= " + children.index + " irank= "
                //+ children.irank + ",idis= " + children.idis + "," + children.cur_list.toString());
        //children.toSolution().write(pw);
        //pw.flush();
        //System.out.println("CrossOX: " + children.cur_list.toString());
        //children.toSolution().print();
        return children;
    }

    //变异 1.Ls,提高才能停止迭代
    //    2.某些点未被加入的突变回来 todo
    public Chromosome mutation(Chromosome chromosome) throws IOException {
        Chromosome new_chromosome = chromosome.copy();
        int iteration = 1;
        //keep still
        //add oneNode
        //delete oneNode
        double time1 = System.currentTimeMillis();
        double[] Probability = new double[2];
        Probability[0] = a1_0 / (a1_0 + a1_1);
        Probability[1] = a1_1 / (a1_0 + a1_1);
        whichlevel = Roulette_wheel_selection(Probability);
        if (whichlevel == 0) {
            new_chromosome = this.changeArtificial2Customer(new_chromosome);
            double time2 = System.currentTimeMillis();
            //System.out.println("A2C cost = "+(time2-time1));
        }
        if (whichlevel == 1) {
            new_chromosome = this.changeCustomer2Artificial(new_chromosome);
            double time2 = System.currentTimeMillis();
            //System.out.println("C2A cost = "+(time2-time1));
        }
        new_chromosome.setFitness1();
        chromosome.setFitness1();
        delta = new_chromosome.fitness < chromosome.fitness;

        Chromosome local_chromosome = this.neighborhood_move(new_chromosome);
        double time3 = System.currentTimeMillis();
        while (iteration < this.UserGAParam.maxLocalCnt) {
            if (local_chromosome.fitness < new_chromosome.fitness) {
                new_chromosome = local_chromosome.copy();
                break;
            } else {

            }
            iteration += 1;
            local_chromosome = this.neighborhood_move(new_chromosome);
        }
        double time4 = System.currentTimeMillis();
        //System.out.println("LocalSearch cost = "+(time4-time3));
        //pw.println("Mutate: new_chromosome_idx = " + new_chromosome.index + " irank= "
                //+ new_chromosome.irank + ",idis= " + new_chromosome.idis);
       // pw.println("level = " + whichlevel + " curlist= "
               // + new_chromosome.cur_list.toString());
        //new_chromosome.toSolution().write(pw);
        //pw.flush();
        updateMutationParam();
        return new_chromosome;
    }

    public void updateMutationParam() {
        if (delta) {
            if (whichlevel == 0) {
                a1_0 += 0.1;
            } else {
                a1_1 += 0.1;
            }
        } else {
            if (whichlevel == 0) {
                a1_0 -= 0.1;
            } else {
                a1_1 -= 0.1;
            }
        }
        a1_0 = Math.max(a1_0, 0);
        a1_1 = Math.max(a1_1, 0);
    }

    public Chromosome changeCustomer2Artificial(Chromosome chromosome) {
        //限制最多不访问的点的个数 把访问的点变成虚拟点

        //List<Integer> visited = chromosome.cur_list.stream().filter(c -> c <= UserConf.N).collect(Collectors.toList());
        List<Integer> artificial = chromosome.cur_list.stream().filter(c -> c > UserConf.N).collect(Collectors.toList());
        if (artificial.size() >= UserGAParam.maxSize) {
            return chromosome;
        }
        int maxC = max(chromosome.cur_list);
        //找一个visited 的点和位置；只要不是N+1之后；找到最大的点和位置；否则不变
        Random r = new Random();
        int i = r.nextInt(this.UserConf.N) + 1;
        int mid = chromosome.cur_list.get(i);
        while (mid >= this.UserConf.N + 1) {
            i = r.nextInt(this.UserConf.N) + 1;
            mid = chromosome.cur_list.get(i);
        }
        //int j = chromosome.cur_list.indexOf(maxC);
        //System.out.println("maxC= " + maxC);
        chromosome.cur_list.set(i, maxC + 1);
        //System.out.println(chromosome.cur_list.toString());
        return chromosome;
        //1- N // N+1 - 2N
    }

    public Chromosome changeArtificial2Customer(Chromosome chromosome) {
        int maxC = max(chromosome.cur_list);
        //if (maxC <= this.UserConf.N) {
        //return chromosome;
        //}
        List<Integer> allC = new ArrayList<>();
        List<Integer> visited = chromosome.cur_list.stream().filter(c -> c <= UserConf.N).collect(Collectors.toList());
        //List<Integer> unvisited = new ArrayList<>();
        for (int k = 1; k <= UserConf.N; k++) {
            allC.add(k);
        }
        //获得没有访问的顾客合集
        allC.removeAll(visited);
        //System.out.println("allC= " + allC.toString());
        Random r = new Random();
        if (allC.size() > 0) {
            int i = r.nextInt(allC.size());
            int j = chromosome.cur_list.indexOf(maxC);
            chromosome.cur_list.set(j, allC.get(i));
        }
        return chromosome;
    }

    //随机交换；两个位置上的点
    public Chromosome neighborhood_move(Chromosome chromosome) throws IOException {
        Random r = new Random();
        Chromosome local = chromosome.copy();
        int i = r.nextInt(this.UserConf.N) + 1;
        int j = r.nextInt(this.UserConf.N) + 1;
        int mid = local.cur_list.get(i);
        local.cur_list.set(i, local.cur_list.get(j));
        local.cur_list.set(j, mid);
        local.setFitness1();
        return local;
    }

    public int Roulette_wheel_selection(double[] Prob) {
        double m = StdRandom.uniform();
        double Probability_Total = 0;
        int ans = 0;
        for (int i = 0; i < Prob.length; i++) {
            Probability_Total += Prob[i];
            if (Probability_Total >= m) {
                ans = i;
                break;
            }
        }
        return ans;
    }

    public void replace() {
        for (int q = 0; q < pop_number; q++) {
            parentPop[q] = tmpPop[q];
            //pw.println("Replace: = " + Arrays.stream(parentPop).map(c -> c.index).collect(Collectors.toList()));
            //pw.flush();
        }
    }

    //根据拥挤度和支配序 排序种群 然后选择前pop_number 个
    public void merge() {
        Chromosome[] newchromosomes = new Chromosome[parentPop.length + tmpPop.length];
        for (int i = 0; i < parentPop.length; i++) {
            newchromosomes[i] = parentPop[i];
        }
        for (int i = 0; i < tmpPop.length; i++) {
            newchromosomes[i + parentPop.length] = tmpPop[i];
        }
        this.quickRank(newchromosomes, false, 0);
        //this.disRank(newchromosomes);
        //pw.write("\n 合并开始前 排序后: ");
        //Arrays.stream(parentPop).forEach(c -> pw.println(c.index + ",idis= " + c.idis + ",irank= " + c.irank));
        //Arrays.stream(tmpPop).forEach(c -> pw.println(c.index + ",idis= " + c.idis + ",irank= " + c.irank));

        //pw.write("\n 合并开始后 排序后: ");
        Collections.sort(Arrays.asList(newchromosomes), MyComparator3);
        for (int q = 0; q < pop_number; q++) {
            parentPop[q] = newchromosomes[q];
        }
        //Arrays.stream(parentPop).map(c -> c.index).collect(Collectors.toList());
        //pw.println("Merge selection: = " + Arrays.stream(parentPop).map(c -> c.index).collect(Collectors.toList()));
        //Arrays.stream(parentPop).forEach(c -> pw.println(c.index + ",idis= " + c.idis + ",irank= " + c.irank));
        //pw.flush();
    }

    // 拥挤度 排序 todo
    public void disRank(Chromosome[] chromosomes) {
        //int[] disN = new int[pop_number];
        //ArrayList<Comparator<Chromosome>> MyComparator = new ArrayList<>();
        Collections.sort(Arrays.asList(chromosomes), MyComparator1);// sort 会改变 Chromosomes的原先的顺序
        Arrays.asList(chromosomes).stream().forEach(c -> c.setFitness1());
        chromosomes[0].idis = UserConf.K;
        chromosomes[chromosomes.length - 1].idis = UserConf.K;
        for (int j = 1; j < chromosomes.length - 1; j++) {
            chromosomes[j].idis += chromosomes[j + 1].chrProfit - chromosomes[j - 1].chrProfit;
        }
        Collections.sort(Arrays.asList(chromosomes), MyComparator2);// sort 是否改变 Chromosomes的原先的顺序
        Arrays.asList(chromosomes).stream().forEach(c -> c.setFitness1());
        chromosomes[0].idis += UserConf.K;
        chromosomes[chromosomes.length - 1].idis += UserConf.K;
        for (int j = 1; j < chromosomes.length - 1; j++) {
            chromosomes[j].idis += chromosomes[j + 1].chrDistance - chromosomes[j - 1].chrDistance;
        }
    }

    // 非支配等级排序
    public void quickRank1(Chromosome[] chromosomes, boolean isWrite, int iter) {
        //计算种群的 np 支配P的个数 SP 被P支配的集合
        for (int i = 0; i < chromosomes.length; i++) {
            //chromosomes[i].index = i;
            chromosomes[i].np = 0;
            chromosomes[i].sp_list = new ArrayList<>();
        }
        Solution[] Solutions = new Solution[chromosomes.length];
        for (int i = 0; i < chromosomes.length; i++) {
            Solutions[i] = chromosomes[i].toSolution();
        }
        //提前先做一遍 不然会很慢
        for (int i = 0; i < chromosomes.length; i++) {
            // for each i
            for (int j = 0; j < chromosomes.length; j++) {
                if (i != j) {
                    Solution Solution_i = Solutions[i];
                    Solution Solution_j = Solutions[j];
                    boolean i2j = Solution_i.isParetoDominate(Solution_j);
                    //System.out.println("i2j = " + i + "," + j + "," + i2j);
                    //boolean j2i = Solution_j.isParetoDominate(Solution_i);
                    //System.out.println("j2i = " + j + "," + i + "," + j2i);
                    if (i2j) {
                        // i 支配 j
                        chromosomes[j].np += 1;
                        chromosomes[i].sp_list.add(chromosomes[j]);
                    }
                /*if (j2i) {
                    // j 支配 i
                    chromosomes[i].np += 1;
                    chromosomes[j].sp_list.add(chromosomes[i]);
                }*/
                }
            }
        }
        List<Chromosome> curF = Arrays.stream(chromosomes).filter(c -> c.np == 0).collect(Collectors.toList());

        if (isWrite) {
            pw2.write("第" + iter + "次迭代的非支配解集\n ");
            curF.stream().forEach(c -> c.toSolution().write(pw2));
            pw2.write("\n ");
            curF.stream().forEach(c -> c.toSolution().write1(pw3, iter));
            pw3.flush();
        }
    }
    public void quickRank(Chromosome[] chromosomes, boolean isWrite, int iter) {
        //计算种群的 np 支配P的个数 SP 被P支配的集合
        for (int i = 0; i < chromosomes.length; i++) {
            //chromosomes[i].index = i;
            chromosomes[i].np = 0;
            chromosomes[i].sp_list = new ArrayList<>();
        }
        Solution[] Solutions = new Solution[chromosomes.length];
        for (int i = 0; i < chromosomes.length; i++) {
            Solutions[i] = chromosomes[i].toSolution();
        }
        //提前先做一遍 不然会很慢
        for (int i = 0; i < chromosomes.length; i++) {
            // for each i
            for (int j = 0; j < chromosomes.length; j++) {
                if (i != j) {
                    Solution Solution_i = Solutions[i];
                    Solution Solution_j = Solutions[j];
                    boolean i2j = Solution_i.isParetoDominate(Solution_j);
                    //System.out.println("i2j = " + i + "," + j + "," + i2j);
                    //boolean j2i = Solution_j.isParetoDominate(Solution_i);
                    //System.out.println("j2i = " + j + "," + i + "," + j2i);
                    if (i2j) {
                        // i 支配 j
                        chromosomes[j].np += 1;
                        chromosomes[i].sp_list.add(chromosomes[j]);
                    }
                /*if (j2i) {
                    // j 支配 i
                    chromosomes[i].np += 1;
                    chromosomes[j].sp_list.add(chromosomes[i]);
                }*/
                }
            }
        }
        List<Chromosome> curF = Arrays.stream(chromosomes).filter(c -> c.np == 0).collect(Collectors.toList());

        if (isWrite) {
            pw2.write("第" + iter + "次迭代的非支配解集\n ");
            curF.stream().forEach(c -> c.toSolution().write(pw2));
            pw2.write("\n ");
            curF.stream().forEach(c -> c.toSolution().write1(pw3, iter));
            pw3.flush();
        }
        curF.forEach(c -> c.irank = 1);
        int sumC = curF.size();
        int rank = 2;

        while (sumC < chromosomes.length) {
            List<Chromosome> curH = new ArrayList<>();
            for (Chromosome c1 : curF) {
                for (Chromosome c2 : c1.sp_list) {
                    c2.np -= 1;
                    if (c2.np == 0) {
                        curH.add(c2);
                    }
                }
            }
            curF = new ArrayList<>();
            curF = curH;
            int finalRank = rank;
            curF.forEach(c -> c.irank = finalRank);
            sumC += curF.size();
            rank += 1;
        }
        System.out.println("quickRank Finish");
    }

    // 遗传算法主流程
    public Chromosome evolve(boolean isNASA) throws IOException {
        System.err.println("Start Initialization ");
        this.bestPop = new Chromosome[this.UserGAParam.maxGeneration];
        double min = 999999;
        double time1 = System.currentTimeMillis();
        this.parentPop = this.initialize();//初始化
        this.quickRank(this.parentPop, true, 0);
        System.err.println("Finish Initialization ");
        for (int i = 0; i < this.UserGAParam.maxGeneration; i++) {
            System.out.println("Generation: = " + i + " start ");
            this.quickRank(this.parentPop, true, i + 1);

            double time2 = System.currentTimeMillis();
           // pw.println("Generation: = " + i + " start ");
            //pw.flush();
            //pw1.println("Generation: = " + i + " start ,每代初始父染色体");
            //Arrays.stream(parentPop).forEach(c -> c.toSolution().write(pw1));
            //pw1.flush();
            System.out.println("Selection start ");
            Chromosome[] mid = this.selection(parentPop);//选择
            double time3 = System.currentTimeMillis();
            System.out.println("Selection finish cost= " + (time3 - time2));
            System.out.println("Crossover start ");
            tmpPop = new Chromosome[pop_number];//子代数组 --- 暂存
            for (int j = 0; j < pop_number; j++) {
                Random r = new Random();
                int fir = r.nextInt(pop_number);
                int sec = r.nextInt(pop_number);
                tmpPop[j] = crossover(mid[fir], mid[sec]);//交叉
            }
            double time4 = System.currentTimeMillis();
            System.out.println("Crossover finish cost= " + (time4 - time3));

            System.out.println("Mutation start ");
            for (int p = 0; p < pop_number; p++) {
                tmpPop[p] = mutation(tmpPop[p]);//变异
            }
            double time5 = System.currentTimeMillis();
            System.out.println("Mutation finish cost= " + (time5 - time4));

            if (isNASA) {
                System.out.println("Merge start ");
                this.merge();
                double time6 = System.currentTimeMillis();
                System.out.println("Merge finish cost= " + (time6 - time5));
            } else {
                System.out.println("Replace start ");
                this.replace();//子代 取代 父类
                double time7 = System.currentTimeMillis();
                System.out.println("Replace finish cost= " + (time7 - time5));
            }
            currentBest = this.getbest(parentPop);
            bestPop[i] = currentBest;
            if (min > currentBest.fitness) {
                bestIndividual = currentBest;
                min = currentBest.fitness;
            }
            System.out.println("CurrentBest Obj= " + min + ":");
            currentBest.toSolution().print();
        }
        // best
        System.out.println("开始收尾");
        Chromosome[] finalChr = new Chromosome[parentPop.length + bestPop.length];
        for (int i = 0; i < parentPop.length; i++) {
            finalChr[i] = parentPop[i];
        }
        for (int i = 0; i < bestPop.length; i++) {
            finalChr[i + parentPop.length] = bestPop[i];
        }
        System.out.println("开始收尾排序 "+finalChr.length);
        this.quickRank1(finalChr, true, this.UserGAParam.maxGeneration + 1);
        finalPop = Arrays.stream(parentPop).filter(c -> c.np == 0).collect(Collectors.toList());
        //最好的和最后一代加起来 存上profit distance 和 路径
        double time8 = System.currentTimeMillis();
        System.out.println("总耗时为 " + (time8 - time1) / 1000.0 + " s");
        return bestIndividual;
    }

    //helper function

    public Chromosome getbest(Chromosome[] chromosomes) {
        Chromosome chromosome = null;
        Arrays.stream(chromosomes).forEach(c -> c.setFitness1());
        //Arrays.stream(chromosomes).forEach(c -> System.out.println("id,obj = " + c.index + "," + c.fitness));
        double min = 99999;
        for (int i = 0; i < pop_number; i++) {
            if (min > chromosomes[i].fitness) {
                min = chromosomes[i].fitness;
                chromosome = chromosomes[i];
            }
        }

        return chromosome;
    }

    public double get_mean(Chromosome[] chromosomes) {
        double ans = 0;
        for (Chromosome chromosome : chromosomes) {
            ans += chromosome.fitness;
        }
        return ans / this.pop_number;
    }

    public Chromosome[] rouletteWheelSelection(Chromosome[] chromosomes) {
        Chromosome[] childrens = new Chromosome[pop_number];

        double[] rouletteWheel;
        //计算适应度
        calRelativeFitness(chromosomes);

        rouletteWheel = new double[this.pop_number];
        rouletteWheel[0] = relativeFitness[0];
        for (int i = 1; i < this.pop_number - 1; i++) {
            rouletteWheel[i] = relativeFitness[i] + rouletteWheel[i - 1];
        }
        rouletteWheel[this.pop_number - 1] = 1;

        for (int i = 0; i < pop_number; i++) {
            double rnd = rand();
            for (int j = 0; j < pop_number; j++) {
                //System.out.println("rnd,prob,j = " + rnd + "," + rouletteWheel[j] + "," + j);
                if (rnd < rouletteWheel[j]) {
                    childrens[i] = chromosomes[j];
                    break;
                }
            }
        }
        return childrens;
    }

    public Chromosome[] eliteselection(Chromosome[] chromosomes) {
        Chromosome[] childrens = new Chromosome[pop_number];
        Random r = new Random();
        for (int t = 0; t < pop_number; t++) {
            int i = r.nextInt(pop_number);
            int j = r.nextInt(pop_number);
            if (chromosomes[i].fitness < chromosomes[j].fitness)
                childrens[t] = chromosomes[i];
            else
                childrens[t] = chromosomes[j];
        }
        return childrens;
    }

    // 计算每个个体的平均适应度
    public double[] calRelativeFitness(Chromosome[] chromosomes) {
        double totalFitness = calTotalFitness(chromosomes);
        for (int i = 0; i < this.pop_number; i++) {
            relativeFitness[i] = Math.abs(chromosomes[i].fitness) / totalFitness;
        }
        //System.out.println("total = " + totalFitness);
        //Arrays.stream(relativeFitness).forEach(c->System.out.println(c));
        return relativeFitness;
    }

    // 计算总体的适应度
    public double calTotalFitness(Chromosome[] chromosomes) {
        double total = 0;
        for (int i = 0; i < chromosomes.length; i++)
            total += Math.abs(chromosomes[i].fitness);
        return total;
    }

    private int rand(int start, int end) {//[start , end)
        return (int) (rand() * (end - start) + start);
    }

    private double rand() {
        return Math.random();
    }

    public static Comparator<Chromosome> MyComparator1 = new Comparator<Chromosome>() {
        @Override
        public int compare(Chromosome e1, Chromosome e2) {
            //升序：前边元素减去后边元素（键升序）
            //降序：后边元素减去前边元素（值降序）
            Solution s1 = e1.toSolution();
            Solution s2 = e2.toSolution();
            s1.getFitness();
            s2.getFitness();
            return (int) (s2.solutionProfit - s1.solutionProfit);
        }
    };
    public static Comparator<Chromosome> MyComparator2 = new Comparator<Chromosome>() {
        @Override
        public int compare(Chromosome e1, Chromosome e2) {
            Solution s1 = e1.toSolution();
            Solution s2 = e2.toSolution();
            s1.getFitness();
            s2.getFitness();
            return (int) (s1.solutionDistance - s2.solutionDistance);
        }
    };
    public static Comparator<Chromosome> MyComparator3 = new Comparator<Chromosome>() {
        @Override
        public int compare(Chromosome e1, Chromosome e2) {
            double nrank1 = e1.irank;
            double nrank2 = e2.irank;
            double nd1 = e1.idis;
            double nd2 = e2.idis;
            int c = 0;
            if (nrank1 < nrank2) {
                c = -1;
            }
            if (nrank1 == nrank2) {
                if (nd1 > nd2) {
                    c = -1;
                } else {
                    if (nd1 == nd2) {
                        c = 0;
                    } else {
                        c = 1;
                    }
                }

            }
            if (nrank1 > nrank2) {
                c = 1;
            }
            return c;
        }
    };

}
