package simulator;

import environment.Cell;
import utils.Constants;
import utils.Utils;

import java.io.*;
import java.util.ArrayList;

/**
 * runnable class from command line
 *
 * @author n.r.zabet@gen.cam.ac.uk
 */
public class SimulatorCLI {

    /**
     * @param args
     * @throws FileNotFoundException
     */
    public static void main(String[] args) throws FileNotFoundException {

        boolean stopAfterBackup = false;
        boolean wasSaved = false;

        Cell cell = null;

        String parametersFilename = "";


        int steps = 10, ensembleSteps;
        double backupAfter;

        if (args.length > 0) {
            parametersFilename = args[0];
        }

        if (args.length > 1) {
            steps = Utils.parseInteger(args[1], 1);
        }

        backupAfter = Constants.MIN_INTERMEDIARY_STATE;

        if (args.length > 2) {
            backupAfter = Utils.parseDouble(args[2], backupAfter);
        }

        if (args.length > 3) {
            stopAfterBackup = Utils.parseBoolean(args[3], stopAfterBackup);
        }


        boolean backupRestarted = false, canRestore;

        int start = 0;
        double time = 0;
        double elapsedTime = 0, estimatedTime;
        int ensemble = 0;
        String lastIntermediaryFilename = "", backupFile, tmpBackupFile;

        //restart a backup file
        File dir = new File(System.getProperty("user.dir"));
        String[] files = dir.list();
        tmpBackupFile = "";
        if (files != null && files.length > 0) {
            int tmp_ensemble = 0, tmp_start = 0, tmp_steps = 0;
            for (String file : files) {
                //located the first backup file
                if (file.startsWith("backup_") && file.endsWith(Constants.INFO_FILE_EXTENSION)) {
                    ArrayList<String> buffer = new ArrayList<String>();
                    String thisLine = "";
                    BufferedReader infoFile = null;
                    canRestore = false;
                    try {
                        infoFile = new BufferedReader(new FileReader(file));
                        while ((thisLine = infoFile.readLine()) != null) {
                            buffer.add(thisLine);
                        }
                        if (buffer.size() >= 4) {
                            //Construct the BufferedWriter object
                            backupFile = buffer.get(0);
                            canRestore = true;
                            if (backupFile != null && !backupFile.isEmpty()) {
                                File f = new File(backupFile);
                                if (!f.exists()) {
                                    f = new File("old_" + backupFile);
                                    if (f.exists()) {
                                        backupFile = "old_" + backupFile;
                                    } else {
                                        canRestore = false;
                                        System.out.println("cannot restore file " + "old_" + backupFile);
                                    }
                                }
                            } else {
                                canRestore = false;
                                System.out.println("cannot restore file " + file);
                            }

                            if (canRestore) {
                                tmp_ensemble = Utils.parseInteger(buffer.get(1), 0);
                                tmp_start = Utils.parseInteger(buffer.get(2), 0);
                                tmp_steps = Utils.parseInteger(buffer.get(3), steps);
                                tmpBackupFile = backupFile;
                            }

                            if (canRestore && (tmp_start) > (start)) {
                                tmpBackupFile = backupFile;
                                lastIntermediaryFilename = backupFile;
                                ensemble = tmp_ensemble;
                                start = tmp_start;
                                steps = tmp_steps;
                            } else {
                                System.out.println("did not restore file " + backupFile);
                            }

                        }
                    } catch (IOException ex) {
                        ex.printStackTrace();
                    }
                }
            }
        }

        if (tmpBackupFile != null && !tmpBackupFile.isEmpty()) {
            FileInputStream fis = null;
            ObjectInputStream in = null;
            try {
                fis = new FileInputStream(tmpBackupFile);
                in = new ObjectInputStream(fis);
                cell = (Cell) in.readObject();
                backupRestarted = true;
                time = cell.cellTime;
                cell.ensemble = ensemble;

                int actualStart = (int) Math.round(cell.cellTime / (cell.totalStopTime / steps));
                System.out.println("start=" + start + "; actualStart = " + actualStart + "");

                if (actualStart != start) {
                    start = actualStart;
                }

                in.close();
                System.out.println("restore file " + tmpBackupFile);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        // FG: parsing of the input parameters file
        if (!backupRestarted) {
            cell = new Cell(parametersFilename, null, true);
        }

        double stopTime = cell.getTotalStopTime();
        double timeStep = stopTime / steps;
        int i = start;

        double lastSave = elapsedTime;

        System.out.println("GRiPE v0.1");

        ensembleSteps = steps * cell.ip.ENSEMBLE_SIZE.value;

        String intermediaryFilename = cell.outputIntermediaryBackupFile;

        //run intervals
        if (cell.totalStopTime > 0) {
            while (i < ensembleSteps && (!stopAfterBackup || !wasSaved)) {
                ensemble = cell.ensemble;
                time += timeStep;
                System.out.println("to perform step " + (i % steps) + " of set " + (i / steps));
                if (i % steps == (steps - 1) || i == ensembleSteps - 1) {
                    time = stopTime;
                    System.out.println("last step of the set to finish");
                }

                elapsedTime += cell.runInterval(time, elapsedTime);
                estimatedTime = steps * elapsedTime / (i + 1) * cell.ip.ENSEMBLE_SIZE.value;
                System.out.println((i + 1) + "/" + ensembleSteps + " elapsed time: " + elapsedTime + "s; estimated " +
						"total time:" + estimatedTime + "s");
                i++;
                if (ensemble != cell.ensemble) {
                    time = 0;
                    System.out.println("set finished");
                }

                //if it has been more than backupAfter seconds from the last save then save the state
                if (backupAfter > 0 && elapsedTime - lastSave > backupAfter) {
                    System.out.println("backup started");
                    double currentTime = System.currentTimeMillis();

                    FileOutputStream fos = null;
                    ObjectOutputStream out = null;
                    try {
                        intermediaryFilename = cell.outputIntermediaryBackupFile;

                        fos = new FileOutputStream("tmp_" + intermediaryFilename);
                        out = new ObjectOutputStream(fos);
                        out.writeObject(cell);
                        out.close();

                        if (!lastIntermediaryFilename.isEmpty()) {
                            File f = new File(lastIntermediaryFilename);
                            if (f.exists() && f.canWrite()) {
                                f.renameTo(new File("old_" + lastIntermediaryFilename));
                            }
                        }

                        // Rename file (or directory)
                        File file = new File("tmp_" + intermediaryFilename);
                        file.renameTo(new File(intermediaryFilename));
                        System.out.println("backup file created " + intermediaryFilename + "; ensemble: " + ensemble + "; step:" + i + "; total:" + steps);

                        BufferedWriter infoFile = null;
                        try {
                            //Construct the BufferedWriter object
                            infoFile = new BufferedWriter(new FileWriter(cell.outputIntermediaryInfoFile));
                            infoFile.write(intermediaryFilename);
                            infoFile.newLine();
                            infoFile.write(ensemble + "");
                            infoFile.newLine();
                            infoFile.write(i + "");
                            infoFile.newLine();
                            infoFile.write(steps + "");
                            infoFile.newLine();
                            infoFile.flush();
                            infoFile.close();
                        } catch (IOException ex) {
                            ex.printStackTrace();
                        }

                        if (!lastIntermediaryFilename.isEmpty()) {
                            File f = new File("old_" + lastIntermediaryFilename);
                            if (f.exists() && f.canWrite()) {
                                f.delete();
                            }
                        }

                        lastIntermediaryFilename = intermediaryFilename;
                        wasSaved = true;
                    } catch (IOException e) {
                        e.printStackTrace();
                    }

                    lastSave = elapsedTime + cell.computeElapsedTime(currentTime);
                }
            }
            //if simulations reached the end the backup files can be deleted
            if (i >= ensembleSteps) {
                if (!lastIntermediaryFilename.isEmpty()) {
                    File f = new File(lastIntermediaryFilename);
                    if (f.exists() && f.canWrite()) {
                        f.delete();
                    }
                }
                File f = new File(cell.outputIntermediaryInfoFile);
                if (f.exists() && f.canWrite()) {
                    f.delete();
                }
            }
        } else {
            cell.runUntilTSReached();
        }
    }

}
