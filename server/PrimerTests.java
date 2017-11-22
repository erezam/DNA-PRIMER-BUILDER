
public class PrimerTests{
    static private int gSum,cSum,aSum,tSum; //


    public static void main(String[] args) {
            String primer = "CCTCTGGAGCGGACTTATTTAC";
            nucleotideSum(primer);
            System.out.println(getGCprecent(primer));
            return;
    }

    private static int getGCprecent(String primer) {   // return G,C precent in primer
        return (int)(((double)(gSum+cSum)/primer.length())*100);
    }
    
  /*  private static int TM(String primer)
    {

    }*/

    private static void nucleotideSum(String primer) // the method counting each nucleotid.
    {
        for(int i=0 ; i<primer.length();i++)
        {
            if(primer.charAt(i)=='A')
                aSum++;
            else if(primer.charAt(i)=='C')
                cSum++;
            else if(primer.charAt(i)=='G')
                gSum++;
            else
                tSum++;
        }
    }

    private static boolean criticalTests(String primer){
        if(primer.charAt(0)=='G')      
            return false;            // primer cant start with G

        int consecutiveG = 0;
        for(int i=0 ; i<primer.length();i++){
            if(primer.charAt(index)=='G')
                consecutiveG++;
            else
                consecutiveG=0;
            if(consecutiveG==3)
                return false;       // primer cant have 3 consecutive G
        }

        

    }
    
}