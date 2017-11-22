
public class PrimerTests{
    static private int gSum,cSum,aSum,tSum; //


    public static void main(String[] args) {
            String primer = "CCTCTGGAGCGGACTTATTTAC";
            System.out.println(GCprecent(primer));
            return;
    }

    private static int getGCprecent(String primer) {   // return G,C precent in primer
        return (int)((double)((gSum+cSum)/primer.length())*100);
    }
    
    private static int TM(String primer)
    {

    }

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
    
}