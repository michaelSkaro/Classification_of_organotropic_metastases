/**
*  Ranks the attributes of a dataset by their gain ratio.
*
* Compilaiton:
*     javac -cp .:/usr/share/java/weka/weka.jar ./src/GainRatio.java -d classes
* Usage:
*     java -cp .:.:/usr/share/java/weka/weka.jar:./classes GainRatio \
*     <Path to dataset (.csv)> <Path of output directory>
*/


import java.io.File;
import java.io.FileWriter;
import java.util.Enumeration;

import weka.core.Instances;
import weka.core.Attribute;
import weka.core.converters.CSVLoader;
import weka.attributeSelection.Ranker;
import weka.attributeSelection.AttributeSelection;
import weka.core.converters.ConverterUtils.DataSource;
import weka.attributeSelection.GainRatioAttributeEval;


public class GainRatio
{
    public static void main(String[] args)
    {
        //for(String arg: args){System.out.println(arg);}
        //System.exit(0);
        String inputPath = args[0];
        String outdirPath = args[1];
        File file = null;
        Instances instances = null;
        CSVLoader loader = new CSVLoader();
        AttributeSelection selector = new AttributeSelection();
        GainRatioAttributeEval eval = new GainRatioAttributeEval();
        Ranker ranker = new Ranker();
        FileWriter out = null;

        try
        {
            file = new File(inputPath);
            loader.setSource(file);
            instances = loader.getDataSet();
            instances.setClassIndex(instances.numAttributes() - 1);

            String[] attributeNames = new String[instances.numAttributes() - 1];
            Enumeration<Attribute> attributeNamesEnum = instances.enumerateAttributes();
            for(int i = 0; i < attributeNames.length; i++)
            {
               attributeNames[i] = attributeNamesEnum.nextElement().name();
            }

            selector.setEvaluator(eval);
            selector.setSearch(ranker);
            selector.SelectAttributes(instances);
            double[][] rankedAttributes = selector.rankedAttributes();

            String outputFilePath = createOutputPath(inputPath, outdirPath);
            out = new FileWriter(outputFilePath);
            out.write("Feature,Score\n");
            for(int i = 0; i < 1000; i++)
            {
                int attributeIndex = (int) rankedAttributes[i][0];
                out.write(attributeNames[attributeIndex] + "," 
                          + rankedAttributes[i][1] + "\n");
            }
            out.close();
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }

    /**
     * Creates the path to the output file.
     *
     * @param inputPath  The file path of the input dataset.
     * @param outdirPath The path of the output directory.
     * @return           The path of the output file.
     */
    private static String createOutputPath(String inputPath, String outdirPath)
    {
        String[] splitPath = inputPath.split("/");
        String fileName = splitPath[splitPath.length - 1];
        String outputFileName = fileName.split("\\.")[0] + "_features_gain_ratio.csv";
        String outputFilePath = null;
        if(outdirPath.charAt(outdirPath.length() - 1) == '/')
        {
            outputFilePath = outdirPath + outputFileName;
        }
        else
        {
            outputFilePath = outdirPath + "/" + outputFileName;
        }
    
        return outputFilePath;
    }
}
