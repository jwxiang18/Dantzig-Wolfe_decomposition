import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class readData {
    public static Instance readdata(String filename) throws IOException {
        // Read the data from the file
        // Return the data as an array of Instance objects
        Instance instance = new Instance();
        try(BufferedReader br = Files.newBufferedReader(Paths.get(filename))){
            instance.NumOrg = Integer.parseInt(br.readLine());
            instance.NumDes = Integer.parseInt(br.readLine());
            instance.NumProduct = Integer.parseInt(br.readLine());


            instance.Supply = new double[instance.NumOrg][instance.NumProduct];
            for (int i = 0; i < instance.NumOrg; i++) {
                String[] line = br.readLine().trim().split("\\s+");
                for (int j = 0; j < instance.NumProduct; j++) {
                    instance.Supply[i][j] = Double.parseDouble(line[j]);
                }
            }

            instance.Demand = new double[instance.NumDes][instance.NumProduct];
            for (int i = 0; i < instance.NumDes; i++) {
                String[] line = br.readLine().trim().split("\\s+");
                for (int j = 0; j < instance.NumProduct; j++) {
                    instance.Demand[i][j] = Double.parseDouble(line[j]);
                }
            }

            instance.Capacity = new double[instance.NumOrg][instance.NumDes];
            for (int i = 0; i < instance.NumOrg; i++) {
                String[] line = br.readLine().trim().split("\\s+");
                for (int j = 0; j < instance.NumDes; j++) {
                    instance.Capacity[i][j] = Double.parseDouble(line[j]);
                }
            }

            instance.Cost = new double[instance.NumOrg][instance.NumDes][instance.NumProduct];
            for (int i = 0; i < instance.NumOrg; i++) {
                for (int j = 0; j < instance.NumDes; j++) {
                    String []line = br.readLine().trim().split("\\s+");
                    for (int k = 0; k < instance.NumProduct; k++) {
                        instance.Cost[i][j][k] = Double.parseDouble(line[k]);
                    }
                }
            }
            return instance;
        }
    }
}
