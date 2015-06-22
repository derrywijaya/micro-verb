package edu.cmu.ml.rtw;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import edu.cmu.ml.rtw.generic.data.annotation.AnnotationType;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.AnnotationTypeNLP;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.AnnotationTypeNLP.Target;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.DependencyParse;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.DependencyParse.Dependency;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.DependencyParse.Node;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.DocumentNLP;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.PoSTag;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.TokenSpan;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.micro.Annotation;
import edu.cmu.ml.rtw.generic.model.annotator.nlp.AnnotatorTokenSpan;
import edu.cmu.ml.rtw.generic.util.Triple;
import edu.cmu.ml.rtw.micro.cat.data.annotation.nlp.AnnotationTypeNLPCat;


public class AnnotationVerb implements AnnotatorTokenSpan<String> {
	private Map<String, Map<String, Double>> verbToRelations, applicableRelations;
	private Map<String, ArrayList<String>> parents, children;
	private Map<String, String> domains;
	private Map<String, String> ranges;

	public static final AnnotationTypeNLP<String> NELL_VERB = new AnnotationTypeNLP<String>("nell-verb", String.class, Target.TOKEN_SPAN);
	private static final AnnotationType<?>[] REQUIRED_ANNOTATIONS = new AnnotationType<?>[] {
		AnnotationTypeNLP.TOKEN,
		AnnotationTypeNLP.SENTENCE,
		AnnotationTypeNLP.POS,
		AnnotationTypeNLP.DEPENDENCY_PARSE,
		AnnotationTypeNLPCat.NELL_CATEGORY,
		AnnotationTypeNLP.LEMMA
	};

	public AnnotationVerb() {
		parents = new HashMap<String, ArrayList<String>>();
		children = new HashMap<String, ArrayList<String>>();
		verbToRelations = new HashMap<String, Map<String, Double>>();
		domains = new HashMap<String, String>();
		ranges = new HashMap<String, String>();
		applicableRelations = new HashMap<String, Map<String, Double>>();
		readHierarchy();
		readMapping();
		readDomainsAndRanges();
		/*PipelineNLPStanford pipelineStanford = new PipelineNLPStanford(30);
		PipelineNLPExtendable pipelineExtendable = new PipelineNLPExtendable();
		pipelineExtendable.extend(new NELLMentionCategorizer());
		PipelineNLP pipeline = pipelineStanford.weld(pipelineExtendable);
		
		DocumentNLP document = new DocumentNLPInMemory(new CatDataTools(), 
				   "Test document", 
				   "Barack Obama was born in Pittsburgh in 2000.",
				   Language.English, pipeline);
		System.out.println("annotating DocumentNLPInMemory");
		annotate(document);*/
		
		/*PipelineNLPStanford pipelineStanford = new PipelineNLPStanford(30);
		System.out.println("Building DocumentNLPInMemory");
		DocumentNLP document = new DocumentNLPInMemory(new DataTools(), 
				   "Test document", 
				   "James was admitted to the hospital in May 2014. He died of cancer in August 2014.",
				   Language.English, pipelineStanford);
		annotate(document);
		*/
		//System.out.println(verbToRelations.get("capture throughout"));
	}

	@Override
	public String getName() {
		return "cmunell_verb-0.0.1";
	}

	@Override
	public boolean measuresConfidence() {
		return true;
	}

	@Override
	public AnnotationType<String> produces() {
		return NELL_VERB;
	}

	@Override
	public AnnotationType<?>[] requires() {
		return REQUIRED_ANNOTATIONS;
	}

	@Override
	public List<Triple<TokenSpan, String, Double>> annotate(DocumentNLP document) {
		Collection<AnnotationTypeNLP<?>> col = new ArrayList<AnnotationTypeNLP<?>>();
		col.add(AnnotationTypeNLPCat.NELL_CATEGORY);
		List<Annotation> annotations = document.toMicroAnnotation(col).getAllAnnotations();	
		Map<Integer, Map<TokenSpan, Map<String, Double>>> billCats = new HashMap<Integer, Map<TokenSpan, Map<String, Double>>>();
		for (Annotation annotation : annotations) {
			//System.out.println(annotation.getSpanStart() + "\t" + annotation.getSpanEnd() + "\t" + annotation.getStringValue());	
			int startTokenIndex = -1;
			int endTokenIndex = -1;
			int sentenceIndex = -1;
			int sentCount = document.getSentenceCount();
			for (int j = 0; j < sentCount; j++) {
				int tokenCount = document.getSentenceTokenCount(j);
				for (int i = 0; i < tokenCount; i++) {
					if (document.getToken(j, i).getCharSpanStart() == annotation.getSpanStart())
						startTokenIndex = i;
					if (document.getToken(j, i).getCharSpanEnd() == annotation.getSpanEnd()) {
						endTokenIndex = i + 1;
						sentenceIndex = j;
						break;
					}
				}
			}
			if (startTokenIndex < 0 || endTokenIndex < 0) {}
			else {
				TokenSpan billCat = new TokenSpan(document, sentenceIndex, startTokenIndex, endTokenIndex); 
				Map<TokenSpan, Map<String, Double>> tempMap = billCats.get(sentenceIndex);
				if (tempMap == null) tempMap = new HashMap<TokenSpan, Map<String, Double>>();
				Map<String, Double> cats = tempMap.get(billCat);
				if (cats == null) cats = new HashMap<String, Double>();
				cats.put(annotation.getStringValue(), annotation.getConfidence());
				tempMap.put(billCat, cats);
				billCats.put(sentenceIndex, tempMap);
			}
		}
		List<Triple<TokenSpan, String, Double>> verbs = new ArrayList<Triple<TokenSpan, String, Double>>();
		
		int sentenceCount = document.getSentenceCount();
		for (int sentIndex = 0; sentIndex < sentenceCount; sentIndex++) {
			DependencyParse sentDP = document.getDependencyParse(sentIndex);
			Node[] sentNodes = sentDP.getTokenNodes();
			for (int nodeIndex = 0; nodeIndex < sentNodes.length; nodeIndex++) {
				Node current = sentNodes[nodeIndex];
				if (current != null) {
					int tokenIndex = current.getTokenIndex();
					PoSTag pos = document.getPoSTag(sentIndex, tokenIndex);
					if (pos.toString().startsWith("V")) {
						List<String> foundVerbs = new ArrayList<String>();
						List<TokenSpan> foundVerbTokenSpans = new ArrayList<TokenSpan>();
						List<Integer> foundVerbObjects = new ArrayList<Integer>();
						
						Map<String, ArrayList<TokenSpan>> foundVerbPPs = new HashMap<String, ArrayList<TokenSpan>>();
						Map<String, ArrayList<Integer>> foundVerbPPObjects = new HashMap<String, ArrayList<Integer>>();
						Map<String, TokenSpan> foundVerbPRTs = new HashMap<String, TokenSpan>();
						ArrayList<Integer> foundDirectObjects = new ArrayList<Integer>();
						ArrayList<Integer> foundSubjects = new ArrayList<Integer>();
						
						boolean isPassive = false;
						String lemmaVerb = document.getTokenAnnotation(AnnotationTypeNLP.LEMMA, sentIndex, tokenIndex);
						TokenSpan tokenVerbSpan = new TokenSpan(document, sentIndex, tokenIndex, tokenIndex + 1);
						Dependency[] dependents = current.getDependents();
						for (Dependency dep : dependents) {
							int depIndex = dep.getDependentTokenIndex();
							String type = dep.getType();
							if (type.startsWith("nmod")) {
								List<Dependency> depsOfNMOD = sentDP.getGovernedDependencies(depIndex);
								for (Dependency depOfNMOD : depsOfNMOD) {
									String typeDepOfNMOD = depOfNMOD.getType();
									if (typeDepOfNMOD.equalsIgnoreCase("case")) {
										int propIndex = depOfNMOD.getDependentTokenIndex();
										String lemmaProp = document.getTokenAnnotation(AnnotationTypeNLP.LEMMA, sentIndex, propIndex);
										String verbPP = lemmaVerb + " " + lemmaProp;
										TokenSpan tokenVerbPPSpan;
										if (tokenIndex < propIndex) {
											tokenVerbPPSpan = new TokenSpan(document, sentIndex, tokenIndex, propIndex + 1);
										} else {
											tokenVerbPPSpan = new TokenSpan(document, sentIndex, propIndex, tokenIndex + 1);
										}
										ArrayList<TokenSpan> tempAL = foundVerbPPs.get(verbPP);
										if (tempAL == null) tempAL = new ArrayList<TokenSpan>();
										tempAL.add(tokenVerbPPSpan);
										ArrayList<Integer> tempALI = foundVerbPPObjects.get(verbPP);
										if (tempALI == null) tempALI = new ArrayList<Integer>();
										tempALI.add(depIndex);
										foundVerbPPs.put(verbPP, tempAL);
										foundVerbPPObjects.put(verbPP, tempALI);
										break;
									}
								}
							} else if (type.startsWith("nsubj") || type.equalsIgnoreCase("agent")) {
								foundSubjects.add(depIndex);
								if (type.equalsIgnoreCase("nsubjpass")) {
									isPassive = true;
								}
							} else if (type.equalsIgnoreCase("prt")) {
								String lemmaPRT = document.getTokenAnnotation(AnnotationTypeNLP.LEMMA, sentIndex, depIndex);
								String verbPRT = lemmaVerb + " " + lemmaPRT;
								TokenSpan tokenVerbPRTSpan;
								if (tokenIndex < depIndex) {
									tokenVerbPRTSpan = new TokenSpan(document, sentIndex, tokenIndex, depIndex + 1);
								} else {
									tokenVerbPRTSpan = new TokenSpan(document, sentIndex, depIndex, tokenIndex + 1);
								}
								foundVerbPRTs.put(verbPRT, tokenVerbPRTSpan);
							} else if (type.equalsIgnoreCase("dobj")) {
								foundDirectObjects.add(depIndex);
							}
						}
						if (foundDirectObjects.size() > 0) {
							if (foundVerbPRTs.size() > 0) {
								for (Map.Entry<String, TokenSpan> e : foundVerbPRTs.entrySet()) {
									String v = e.getKey();
									TokenSpan t = e.getValue();
									for (int obj : foundDirectObjects) {
										foundVerbs.add(v);
										foundVerbTokenSpans.add(t);
										foundVerbObjects.add(obj);		
									}
								}
							} else {
								for (int obj : foundDirectObjects) {
									foundVerbs.add(lemmaVerb);
									foundVerbTokenSpans.add(tokenVerbSpan);
									foundVerbObjects.add(obj);									
								}
							}
						}
						if (foundVerbPPs.size() > 0) {
							for (Map.Entry<String, ArrayList<TokenSpan>> e : foundVerbPPs.entrySet()) {
								String v = e.getKey();
								ArrayList<TokenSpan> ts = e.getValue();
								ArrayList<Integer> objs = foundVerbPPObjects.get(v);
								for (int ii = 0; ii < ts.size(); ii++) {
									TokenSpan t = ts.get(ii);
									Integer obj = objs.get(ii);
									foundVerbs.add(v);
									foundVerbTokenSpans.add(t);
									foundVerbObjects.add(obj);
								}
							}
						}
						for (int i = 0; i < foundVerbs.size(); i++) {
							for (int subj : foundSubjects) {
								int obj = foundVerbObjects.get(i);
								String lemmaV = foundVerbs.get(i);
								TokenSpan tokenSpanV = foundVerbTokenSpans.get(i);
								Map<TokenSpan, Map<String, Double>> tempMap = billCats.get(sentIndex);
								if (isPassive) lemmaV = "(passive) " + lemmaV.trim();
								ArrayList<String> catSubj = null; ArrayList<String> catObj = null;
								ArrayList<String> topCatSubj = null; ArrayList<String> topCatObj = null;
								if (tempMap != null) {
									for (Map.Entry<TokenSpan, Map<String, Double>> e : tempMap.entrySet()) {
										TokenSpan ct = e.getKey();
										if (ct.containsToken(sentIndex, subj)) {
											catSubj = new ArrayList<String>();
											double maxVal = -1;
											for (Map.Entry<String, Double> ent : e.getValue().entrySet()) {
												catSubj.add(ent.getKey());
												double val = ent.getValue();
												if (val > maxVal) {
													if (topCatSubj == null) topCatSubj = new ArrayList<String>();
													maxVal = val;
													topCatSubj.clear();
													topCatSubj.add(ent.getKey());
												} else if (val == maxVal) {
													topCatSubj.add(ent.getKey());
												}
											}
										}
										if (ct.containsToken(sentIndex, obj)) {
											catObj = new ArrayList<String>();
											double maxVal = -1;
											for (Map.Entry<String, Double> ent : e.getValue().entrySet()) {
												catObj.add(ent.getKey());
												double val = ent.getValue();
												if (val > maxVal) {
													if (topCatObj == null) topCatObj = new ArrayList<String>();
													maxVal = val;
													topCatObj.clear();
													topCatObj.add(ent.getKey());
												} else if (val == maxVal) {
													topCatObj.add(ent.getKey());
												}
											}
										}
									}
									
								}
								//System.out.println(catSubj.toString() + "\t" + catObj.toString() + "\t" + lemmaV);
								if (catSubj != null && catObj != null) {
									if (verbToRelations.get(lemmaV) != null) {
										Map<String, Double> candidateRelations = verbToRelations.get(lemmaV);
										for (Map.Entry<String, Double> entry : candidateRelations.entrySet()) {
											String domain = entry.getKey().split(" ")[0];
											String relation = entry.getKey().split(" ")[1];
											String range = entry.getKey().split(" ")[2];
											double conf = entry.getValue();
											if (hasCategory(catSubj, domain) && hasCategory(catObj, range)) {
												//System.out.println(lemmaSubj + "(" + catSubj + ")" + "\t" + lemmaV + "\t" + lemmaObj + "(" + catObj + ")");
												verbs.add(new Triple<TokenSpan, String, Double>
													(new TokenSpan(document, sentIndex, tokenSpanV.getStartTokenIndex(), tokenSpanV.getEndTokenIndex()), 
															relation, conf));
											}
										}
									} else {
										for (String s : topCatSubj) {
											for (String o : topCatObj) {
												String type = s + ":::" + o;
												Map<String, Double> candidateRels = applicableRelations.get(type);
												if (candidateRels != null) {
													Iterator<Map.Entry<String, Double>> it1 = candidateRels.entrySet().iterator();
													while (it1.hasNext()) {
														Map.Entry<String, Double> ent1 = it1.next();
														String relation = ent1.getKey();
														double conf = ent1.getValue();
														verbs.add(new Triple<TokenSpan, String, Double>
														(new TokenSpan(document, sentIndex, tokenSpanV.getStartTokenIndex(), tokenSpanV.getEndTokenIndex()), 
																relation, conf));														
													}
												}
											}
										}
										
									}
								} else {
									if (verbToRelations.get(lemmaV) != null) {
										Map<String, Double> candidateRelations = verbToRelations.get(lemmaV);
										for (Map.Entry<String, Double> entry : candidateRelations.entrySet()) {
											String relation = entry.getKey().split(" ")[1];
											double conf = entry.getValue();
											verbs.add(new Triple<TokenSpan, String, Double>
											(new TokenSpan(document, sentIndex, tokenSpanV.getStartTokenIndex(), tokenSpanV.getEndTokenIndex()), 
													relation, conf));
										}
									}
								}
							}
						}
					}
				}
			}
		}
		return verbs;
	}

	private boolean hasCategory(ArrayList<String> catSubj, String domain) {
		ArrayList<String> parent = parents.get(domain);
		if (catSubj.contains(domain)) return true;
		for (String p : parent) {
			if (catSubj.contains(p)) return true;
		}
		return false;
	}

	private void readHierarchy() {
		try {
			InputStream is = AnnotationVerb.class.getResourceAsStream("/NELL-hierarchy.txt");
			BufferedReader bfr = new BufferedReader(new InputStreamReader(is));
			String line, temp[];
			while ((line = bfr.readLine()) != null) {
				temp = line.split("\t");
				ArrayList<String> tempAL = parents.get(temp[0].trim().replace("concept:", ""));
				if (tempAL == null) tempAL = new ArrayList<String>();
				tempAL.add(temp[2].trim().replace("concept:", ""));
				parents.put(temp[0].trim().replace("concept:", ""), tempAL);
				
				ArrayList<String> tempCL = children.get(temp[2].trim().replace("concept:", ""));
				if (tempCL == null) tempCL = new ArrayList<String>();
				tempCL.add(temp[0].trim().replace("concept:", ""));
				children.put(temp[2].trim().replace("concept:", ""), tempCL);
			}
			bfr.close();
			
		} catch (IOException ioe) {
		    throw new RuntimeException(ioe);
		}
	}
		
	private void readDomainsAndRanges() {
		try {
			InputStream is = AnnotationVerb.class.getResourceAsStream("/NELL-relationTypes.txt");
			BufferedReader bfr = new BufferedReader(new InputStreamReader(is));
			String line, temp[];
			while ((line = bfr.readLine()) != null) {
				temp = line.toLowerCase().split("\t");
				String relation = temp[0].trim().replace("concept:", "");
				String kind = temp[1].trim();
				String type = temp[2].trim().replace("concept:", "");
				if (kind.equalsIgnoreCase("domain")) {
					domains.put(relation, type);
				} else if (kind.equalsIgnoreCase("range")) {
					ranges.put(relation, type);
				}
			}
			bfr.close();
		} catch (IOException ioe) {
		    throw new RuntimeException(ioe);
		}

		Map<String, Double> relationConfs = new HashMap<String, Double>();
		try {
			InputStream is = AnnotationVerb.class.getResourceAsStream("/priorProbRelation.txt");
			BufferedReader bfr = new BufferedReader(new InputStreamReader(is));
			String line, temp[];
			while ((line = bfr.readLine()) != null) {
				temp = line.split("\t");
				String relation = temp[1].trim().replace("concept:", "");
				double conf = Double.parseDouble(temp[2].trim());
				relationConfs.put(relation, conf);
			}
			bfr.close();
		} catch (IOException ioe) {
		    throw new RuntimeException(ioe);
		}

		Iterator<Map.Entry<String, String>> it = domains.entrySet().iterator();
		while (it.hasNext()) {
			Map.Entry<String, String> e = it.next();
			String relation = e.getKey();
			if (relationConfs.get(relation) != null) {
				double conf = relationConfs.get(relation);
				String domain = e.getValue();
				String range = ranges.get(relation);
				
				ArrayList<String> childrenDomain = children.get(domain);
				ArrayList<String> childrenRange = children.get(range);
				if (childrenDomain == null) childrenDomain = new ArrayList<String>();
				if (childrenRange == null) childrenRange = new ArrayList<String>();
				childrenDomain.add(domain);
				childrenRange.add(range);
				
				for (String d : childrenDomain) {
					for (String r : childrenRange) {
						String type = d + ":::" + r;
						Map<String, Double> tempAL = applicableRelations.get(type);
						if (tempAL == null) tempAL = new HashMap<String, Double>();
						if (tempAL.size() == 0) tempAL.put(relation, conf);
						else {
							Iterator<Map.Entry<String, Double>> it1 = tempAL.entrySet().iterator();
							double prevVal = it1.next().getValue();
							if (conf > prevVal) {
								tempAL.clear();
								tempAL.put(relation, conf);
							} else if (conf == prevVal) {
								tempAL.put(relation, conf);
							}
						}
						applicableRelations.put(type, tempAL);
					}
				}				
			}
		}
	}
	
	private void readMapping() {
		try {
			InputStream is = AnnotationVerb.class.getResourceAsStream("/naive-bayes-em-mapping-filtered-with-type-checking-and-prior-and-cps-with-decay-3-knee-verbs.txt");
			BufferedReader bfr = new BufferedReader(new InputStreamReader(is));
			String line;
			String temp[], relations[], verbs[], domain, range, relation, verb;
			double score;
			
			while ((line = bfr.readLine()) != null) {
				temp = line.split("\t");
				relations = temp[0].split(" ");
				domain = relations[0].trim().replace("concept:", "");
				range = relations[2].trim().replace("concept:", "");
				relation = relations[1].trim().replace("concept:", "");
				temp = temp[1].split(";");
				for (String t : temp) {
					verbs = t.split(",");
					verb = verbs[0].trim();
					verb = verb.replace("(also in CPL)", "").trim();
					score = Double.parseDouble(verbs[1].trim());
					String rrelation;
					if (verb.startsWith("arg1")) {
						rrelation = domain + " " + relation + " " + range;
						verb = verb.replace("arg1 ", "").trim();
						verb = verb.replace(" arg2", "").trim();
					} else {
						rrelation = range + " " + relation + " " + domain;
						verb = verb.replace("arg2 ", "").trim();
						verb = verb.replace(" arg1", "").trim();
					}
					Map<String, Double> tempMap = verbToRelations.get(verb);
					if (tempMap == null) tempMap = new HashMap<String, Double>();
					tempMap.put(rrelation, score);
					verbToRelations.put(verb, tempMap);
				}
			}
			bfr.close();
			
		} catch (IOException ioe) {
		    throw new RuntimeException(ioe);
		}
	}
}