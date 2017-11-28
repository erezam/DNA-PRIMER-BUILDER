
public class PrimerTests{

    public static void main(String[] args) {
            String forwardPrimer = "CTACCTCCACCATGCCAAGT";
            String reversePrimer = "CCTCTGGAGCGGACTTATTTAC";
            Primer p1 = new Primer(forwardPrimer);
            Primer p2 = new Primer(reversePrimer);
            System.out.println("Forward tm:"+p1.getTm());
            System.out.println("Reverse tm:"+p2.getTm());
            System.out.println("Forward GC:"+p1.getPercentCg());
            System.out.println("Reverse GC:"+p2.getPercentCg());
            /*
            System.out.println("forwardPrimer syntax: "+syntaxTests(p1.getPrimerStr()));
            System.out.println("reversePrimer syntax: "+syntaxTests(p2.getPrimerStr()));
            System.out.println("forwardPrimer tmTest: "+ tmTest(p1));
            System.out.println("reversePrimer tmTest: "+ tmTest(p2));
            System.out.println("temperture Differences: " + tmDifferencesTest(p1,p2));*/
            return;
    }

    private static boolean syntaxTests(String primer){
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

    private static boolean tmTest(Primer primer)
    {
        if(primer.getTm()< 50 ||primer.getTm() > 70 )
             return false;
        return true;
    }
    
    private static boolean tmDifferencesTest(Primer primerF , Primer primerR)
    {
        if(Math.abs(primerF.getTm() - primerR.getTm())>5) 
             return false;
        return true;
    }
}