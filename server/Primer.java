class Primer{
    private int gSum,cSum,aSum,tSum,tm,prcentGC,Length;
    private String primerStr;

    public Primer(String primer){
        this.primerStr = primer;
        for(int i=0 ; i<primer.length();i++)
        {
            if(primer.charAt(i)=='A')
                this.aSum++;
            else if(primer.charAt(i)=='C')
                this.cSum++;
            else if(primer.charAt(i)=='G')
                this.gSum++;
            else
                this.tSum++;
        }
        this.prcentGC=(int)(((double)(gSum+cSum)/primer.length())*100);
        this.tm=(int)(64.9 + 41*(gSum+cSum-16.4)/(aSum+tSum+gSum+cSum));
        this.Length = primer.length();
    }

    public String getPrimerStr(){
        return this.primerStr;
    }
    
    public int getgSum(){
        return this.gSum;
    }

    public int getaSum(){
        return this.aSum;
    }

    public int gettSum(){
        return this.tSum;
    }

    public int getcSum(){
        return this.cSum;
    }

    public int getTm(){
        return this.tm;
    }

    public int getPercentCg(){
        return this.prcentGC;
    }

    public int getLength(){
        return this.Length;
    }
   
}