/**
 * Constructs a Random Forest, performs 10-Fold 
 * cross-validation, and saves the results.
 * Intended for pairwise classifications.
 *
 * Usage:
 * java -cp .:./lib/*:classes CancerClassifiers <input_dir>/ <output_dir>/
 *
 * Note: All supplied directory paths must end with '/', 
 * or else an error will be thrown.
 */

import java.io.File;
import java.lang.Double;
import java.util.Random;
import java.nio.file.Path;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;

import weka.core.Instances;
import weka.classifiers.trees.J48;
import weka.classifiers.Evaluation;
import weka.classifiers.trees.RandomForest;
import weka.core.converters.ConverterUtils.DataSource;

public class CancerClassifiers
{
    public static void main(String[] args) 
    {
        try
        {
            // Read in an input directory from the command line
            String inputDirString = args[0]; 
            File inputDir = new File(inputDirString);
            String[] fileNames = inputDir.list();

            //Read in an output directory from the command line
            String outputDirString = args[1];
            File outputDir = new File(outputDirString);
            Path outputDirPath = outputDir.toPath();
            if(!Files.isDirectory(outputDirPath))
            {
                Files.createDirectory(outputDirPath);
            }

            File outputFile = null;
            Path outputPath = null;
            for(String fileName: fileNames)
            {
                System.out.println(fileName);
                String[] names = fileName.split("TCGA-");
                //String output_name = names[1] + names[2].split("\\.")[0];
                String output_name = fileName;
                
                //Create Instances
                DataSource source = new DataSource(inputDirString + fileName);
                Instances instances = source.getDataSet();
                instances.setClassIndex(instances.numAttributes() - 1);
            
                //10-Fold cross-validation evaluations
                RandomForest tree = new RandomForest();
                Evaluation eval = new Evaluation(instances);
                eval.crossValidateModel(tree, instances, 10, new Random());
                String summary = eval.toSummaryString();
                String classDetails = eval.toClassDetailsString();
                String confusionMatrix = eval.toMatrixString();

                //Create and write to an output file
                String outputString = outputDir + "/" + output_name + ".txt";
                outputFile = new File(outputString);
                outputPath = outputFile.toPath();
                Files.writeString(outputPath, summary + '\n'
                                    + classDetails + '\n'
                                    + confusionMatrix);
            }
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }
    
}
