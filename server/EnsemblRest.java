package transcript;
import java.net.URL;
import java.net.URLConnection;
import java.net.HttpURLConnection;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.io.Reader;
 
 
public class EnsemblRest {
 
  public static void main(String[] args) throws Exception {
      String id = getGeneTranscriptId("APOE");
      //String cDna = getCdna(id);
      System.out.println("tran id: "+id);
      //System.out.println(cDna);
  }

  private static String getGeneTranscriptId(String geneName) throws Exception {
    String server = "https://rest.ensembl.org";
    String ext = "/lookup/symbol/homo_sapiens/"+geneName+"?expand=1";
    URL url = new URL(server + ext);
 
    URLConnection connection = url.openConnection();
    HttpURLConnection httpConnection = (HttpURLConnection)connection;
    
    httpConnection.setRequestProperty("Content-Type", "application/json");
 
    InputStream response = connection.getInputStream();
    int responseCode = httpConnection.getResponseCode();
 
    if(responseCode != 200) {
      throw new RuntimeException("Response code was not 200. Detected response was "+responseCode);
    }
 
    String output;
    Reader reader = null;
    try {
      reader = new BufferedReader(new InputStreamReader(response, "UTF-8"));
      StringBuilder builder = new StringBuilder();
      char[] buffer = new char[8192];
      int read;
      while ((read = reader.read(buffer, 0, buffer.length)) > 0) {
        builder.append(buffer, 0, read);
      }
      output = builder.toString();
    } 
    finally {
        if (reader != null) try {
          reader.close(); 
        } catch (IOException logOrIgnore) {
          logOrIgnore.printStackTrace();
        }
    }
    System.out.println(output);
    int index = output.indexOf("id", output.indexOf(geneName+"-201"))+5;
    String transcriptId= output.substring(index, index+15);
    return transcriptId;
  }

  private static String getCdna(String geneID) throws Exception {
    String server = "http://rest.ensembl.org";
    String ext = "/sequence/id/ENST00000288602?type=cdna";
    URL url = new URL(server + ext);
 
    URLConnection connection = url.openConnection();
    HttpURLConnection httpConnection = (HttpURLConnection)connection;
    
    httpConnection.setRequestProperty("Content-Type", "text/x-fasta");
    
 
    InputStream response = connection.getInputStream();
    int responseCode = httpConnection.getResponseCode();
 
    if(responseCode != 200) {
      throw new RuntimeException("Response code was not 200. Detected response was "+responseCode);
    }
 
    String output;
    Reader reader = null;
    try {
      reader = new BufferedReader(new InputStreamReader(response, "UTF-8"));
      StringBuilder builder = new StringBuilder();
      char[] buffer = new char[8192];
      int read;
      while ((read = reader.read(buffer, 0, buffer.length)) > 0) {
        builder.append(buffer, 0, read);
      }
      output = builder.toString();
    } 
    finally {
        if (reader != null) try {
          reader.close(); 
        } catch (IOException logOrIgnore) {
          logOrIgnore.printStackTrace();
        }
    }
 
    return output;
  }
}