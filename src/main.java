import gurobi.GRBException;

import java.io.IOException;

public class main {
    public static void main(String[] args) throws IOException, GRBException {
        Instance instance = readData.readdata("multicommodity.dat");
        DW dw = new DW(instance);
        dw.main();
    }
}
