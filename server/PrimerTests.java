
public class PrimerTests{
    static private int gSum,cSum,aSum,tSum; //


    public static void main(String[] args) {
            String primer = "CCTCTGGAGCGGACTTATTTAC";
            nucleotideSum(primer);
            System.out.println(criticalTests(primer));
            System.out.println(getGCprecent(primer));
            System.out.println(getTM(primer));
            return;
    }

    private static int getGCprecent(String primer) {   // return G,C precent in primer
        return (int)(((double)(gSum+cSum)/primer.length())*100);
    }
    
    private static int getTM(String primer)
    {
            int temperture = (int)(64.9 + 41*(gSum+cSum-16.4)/(aSum+tSum+gSum+cSum));
            return temperture;
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

    private static boolean criticalTests(String primer){
        if(primer.charAt(0)=='G') 
        {
            System.out.println("primer cant start with G");
            return false;            // primer cant start with G
        }
        if(primer.contains("GGG"))
        {
            System.out.println("primer cant have 3 consecutive G");
            return false;        // primer cant have 3 consecutive G
        }
        if(primer.contains("GGAG"))
        {
            System.out.println("primer cant have GGAG");
            return false;        // primer cant have GGAG
        }
        if(primer.contains("AAAAAA"))
        {
            System.out.println("primer cant have AAAAAA");
            return false;       // primer cant have 6 consecutive A
        }   
        int consecutiveC = 0;
        for(int i=1 ; i<primer.length()-1;i++){
            if(primer.charAt(i)=='C')
                consecutiveC++;
            else
                consecutiveC=0;
            if(consecutiveC==2)
            {
                System.out.println("primer cant have cc in the middle");
                return false;       // primer cant have 2 consecutive C in the middle                
            }
        }

        return true;

        

    }
    
}