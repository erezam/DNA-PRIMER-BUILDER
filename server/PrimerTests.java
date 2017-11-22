
public class PrimerTests{

    public static void main(String[] args) {
            
            String forwardPrimer = "CCTCTGGGGAGCGGACTTATTTAC";
            String reversePrimer = "CCTCTGGAGCGGACTTATTTAC";
            Primer p1 = new Primer(forwardPrimer);
            System.out.println(criticalTests(p1.getPrimerStr()));
            return;
    }

    private static boolean criticalTests(String primer){
        if(primer.charAt(0)=='G') 
        {
            System.out.println("primer cant start with G");
            return false;            // primer cant start with G
        }
        if(primer.contains("GGGG"))
        {
            System.out.println("primer cant have 4 consecutive G");
            return false;        // primer cant have 3 consecutive G
        }
        if(primer.contains("CCCC"))
        {
            System.out.println("primer cant have 4 consecutive C");
            return false;        // primer cant have 4 consecutive C
        }
        if(primer.contains("AAAAAA"))
        {
            System.out.println("primer cant have AAAAAA");
            return false;       // primer cant have 6 consecutive A
        }   
        /*int consecutiveC = 0;
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
        }*/
        return true;
    }
    
}