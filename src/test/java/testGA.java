import Algo.Chromosome;
import Algo.GA_strategy;
import base.Customer;
import org.junit.Test;
import param.Conf;
import param.GAParam;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

public class testGA {

    public static void main(String[] args) throws IOException {
        FileWriter fw1 = new FileWriter("./dataset/probLog.txt");
        BufferedWriter bw1 = new BufferedWriter(fw1);
        PrintWriter pw1 = new PrintWriter(bw1);
        //里程的设置
        GAParam UserGAParam = new GAParam();
        Conf UserCof = new Conf();
        int a = 70, c = 90, b = 110;
        int CustomerNumber = 25;
        boolean isNASA = true;
        for (int k = 1; k < 9; k++) {//RC108-208 R101-R112 R201-R211 C208 C109
            String prefix = "C" + Integer.toString(200 + k) + "_25";
            System.err.println("Solve "+prefix);
            UserCof.readInstance("./dataset/inputdata/" + prefix + ".txt", a, c, b, pw1);
            GA_strategy GA = new GA_strategy(UserCof, UserGAParam, prefix);

            System.out.println("运行中 ");
            Chromosome best = GA.evolve(isNASA);
            System.out.println("最优解 ");
            best.toSolution().print();

            FileWriter fw = new FileWriter("./dataset/" + prefix + "result.txt");
            BufferedWriter bw = new BufferedWriter(fw);
            PrintWriter pw = new PrintWriter(bw);

            pw.println("\n单目标 Best: chromosome_idx = " + best.index + " irank= "
                    + best.irank + ",idis= " + best.idis);
            best.toSolution().write(pw);
            pw.flush();

            pw.print("\n多目标种群 Current_best_pop:");
            for (int i = 0; i < GA.bestPop.length; i++) {
                pw.print("\nBest pop: " + i + " : " + GA.bestPop[i].index);
                GA.bestPop[i].toSolution().write(pw);
                pw.flush();
            }
            pw.println("");

            pw.print("\n最后一代种群混合最好的 :");
            GA.finalPop.stream().forEach(p -> p.toSolution().write(pw));
            pw.flush();
            pw.println("");
            pw.close();
        }

    }
}
