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
    String server = "https://rest.ensembl.org";
    String ext = "/sequence/id/ENSCAFT00000046249.1?type=cdna";
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
 
    System.out.println(output);
  }
}