
public class PrimerTests{
    static private int gSum,cSum,aSum,tSum; //


    public static void main(String[] args) {
            String primer = "CCTCTGGAGCGGACTTATTTAC";
            System.out.println(GCprecent(primer));
            return;
    }

    private static int GCprecent(String primer) {   // return G,C precent in primer
        double count=0;
        for(int i=0 ; i<primer.length();i++)
            if((primer.charAt(i)=='G') || (primer.charAt(i)=='C'))
                count+=1; 
        System.out.println(primer.length());        
        int tmp = (int)((count/(double)primer.length())*100);
        return tmp;
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