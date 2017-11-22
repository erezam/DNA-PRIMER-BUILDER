import java.util.Scanner;

public class PrimerTests{
    public static void main(String[] args) {
            Scanner input = new Scanner(System.in);
            System.out.println("Enter Primer:");
            String primer;
            primer = input.next();
            //System.out.println(GCprecent(primer));
            return;
    }

    private static int GCprecent(String primer) {   // return G,C precent in primer
        int count=0;
        for(int i=0 ; i<primer.length();i++)
            if((primer.charAt(i)=='G') || (primer.charAt(i)=='C'))
                count++; 
        return (count/primer.length())*100;
    }
}