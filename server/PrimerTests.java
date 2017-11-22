
public class PrimerTests{
    public static void main(String[] args) {
            String primer = "CCTCTGGAGCGGACTTATTTAC";
            System.out.println(GCprecent(primer));
            return;
    }

    private static int getGCprecent(String primer) {   // return G,C precent in primer
        return (int)((double)((gSum+cSum)/primer.length())*100);
    }

    
}