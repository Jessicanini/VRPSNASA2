package util;

//helper function
public class utils {

    public utils() {
    }

    public static boolean greater(double a, double b) {
        if (a - b > 1E-6) {
            return true;
        } else {
            return false;
        }
    }

    public static boolean greatereq(double a, double b) {
        if (a - b >= 1E-6 || a==b ) {
            return true;
        } else {
            return false;
        }
    }


}
