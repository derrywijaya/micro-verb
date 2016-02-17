package edu.cmu.ml.rtw;

import java.util.List;

import edu.cmu.ml.rtw.generic.data.DataTools;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.DocumentNLP;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.DocumentNLPInMemory;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.DocumentNLPMutable;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.Language;
//import edu.cmu.ml.rtw.generic.data.annotation.nlp.TokenSpan;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.micro.Annotation;
import edu.cmu.ml.rtw.generic.model.annotator.nlp.PipelineNLP;
import edu.cmu.ml.rtw.generic.model.annotator.nlp.PipelineNLPExtendable;
import edu.cmu.ml.rtw.generic.model.annotator.nlp.PipelineNLPStanford;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.SerializerDocumentNLPMicro;
//import edu.cmu.ml.rtw.generic.util.Triple;
//import edu.cmu.ml.rtw.micro.cat.data.CatDataTools;
import edu.cmu.ml.rtw.micro.cat.data.annotation.nlp.NELLMentionCategorizer;
import edu.cmu.ml.rtw.contextless.ContextlessNPCategorizer;

/**
 * Hello world!
 *
 */
public class App 
{
    public static void main( String[] args )
    {
        System.out.println( "Testing Verb Annotation!" );
        PipelineNLPStanford pipelineStanford = new PipelineNLPStanford(30);
        PipelineNLPExtendable pipelineExtendable = new PipelineNLPExtendable();
        pipelineExtendable.extend(new NELLMentionCategorizer());
        pipelineExtendable.extend(new ContextlessNPCategorizer());
        pipelineExtendable.extend(new AnnotationVerb());
        PipelineNLP pipeline = pipelineStanford.weld(pipelineExtendable);
        DataTools dataTools = new DataTools();
        dataTools.addAnnotationTypeNLP(AnnotationVerb.NELL_VERB);
        DocumentNLPMutable document = new DocumentNLPInMemory(dataTools, 
          "Test document", 
          "Barack Obama called his wife. Barack Obama called Michelle Obama.");
        pipeline.run(document);
        SerializerDocumentNLPMicro microSerial = new SerializerDocumentNLPMicro(dataTools);
        List<Annotation> annotations = microSerial.serialize(document).getAllAnnotations();
        for (Annotation annotation : annotations) {
          System.out.println(annotation.toJsonString());
        }
        
        
        /*PipelineNLPStanford pipelineStanford = new PipelineNLPStanford(30);
		PipelineNLPExtendable pipelineExtendable = new PipelineNLPExtendable();
		pipelineExtendable.extend(new NELLMentionCategorizer());
		PipelineNLP pipeline = pipelineStanford.weld(pipelineExtendable);
		
		DocumentNLP document = new DocumentNLPInMemory(new CatDataTools(), 
				   "Test document", 
				   "Barack Obama was born in Pittsburgh in 2000.");
		System.out.println("annotating DocumentNLPInMemory");
                pipeline.run(document);
		List<Triple<TokenSpan, String, Double>> verbs = av.annotate(document);
		for (Triple<TokenSpan, String, Double> e : verbs) {
			System.out.println("sentenceIndex = " + e.getFirst().getSentenceIndex() + "\t" + 
					"startToken = " + e.getFirst().getStartTokenIndex() + "\t" +
					"endToken = " + e.getFirst().getEndTokenIndex() + "\t" +
					"verbAnnotation = " + e.getSecond() + "\t" + "verbConfidence = " + e.getThird());
		}*/
    }
}
