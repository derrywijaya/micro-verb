package edu.cmu.ml.rtw;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
//import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import edu.cmu.ml.rtw.generic.data.annotation.AnnotationType;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.AnnotationTypeNLP;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.AnnotationTypeNLP.Target;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.DependencyParse;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.DependencyParse.Dependency;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.DependencyParse.Node;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.DocumentNLP;
//import edu.cmu.ml.rtw.generic.data.annotation.nlp.DocumentNLPMutable;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.PoSTag;
import edu.cmu.ml.rtw.generic.data.annotation.nlp.TokenSpan;
//import edu.cmu.ml.rtw.generic.data.annotation.nlp.micro.Annotation;
import edu.cmu.ml.rtw.generic.model.annotator.nlp.AnnotatorTokenSpan;
import edu.cmu.ml.rtw.generic.util.Triple;
import edu.cmu.ml.rtw.micro.cat.data.annotation.nlp.AnnotationTypeNLPCat;
//import edu.cmu.ml.rtw.micro.cat.data.CatDataTools;
import edu.cmu.ml.rtw.contextless.ContextlessNPCategorizer;
import edu.cmu.ml.rtw.contextless.TypedNP;


public class AnnotationVerb implements AnnotatorTokenSpan<String> {
	private Map<String, Map<String, Double>> verbToRelations;//, applicableRelations;
	private Map<String, Map<String, Double>> typeToRelations;
	private Map<String, Map<String, Integer>> parents;
	private Map<String, Double> relationPriors;
	//private int LIMIT = 3; // limit number of relations returned
	//private Map<String, ArrayList<String>> parents, children;
	//private Map<String, Vector<String>> parentsSorted;
	//private Map<String, ArrayList<String>> children;
	//private Map<String, String> domains;
	//private Map<String, String> ranges;

	public static final AnnotationTypeNLP<String> NELL_VERB = new AnnotationTypeNLP<String>("nell-verb", String.class, Target.TOKEN_SPAN);
	private static final AnnotationType<?>[] REQUIRED_ANNOTATIONS = new AnnotationType<?>[] {
		AnnotationTypeNLP.TOKEN,
		AnnotationTypeNLP.SENTENCE,
		AnnotationTypeNLP.POS,
		AnnotationTypeNLP.DEPENDENCY_PARSE,
		AnnotationTypeNLPCat.NELL_CATEGORY,
		AnnotationTypeNLP.LEMMA,
		ContextlessNPCategorizer.OUTOFCONTEXT_NP_CATEGORIES
	};

	public AnnotationVerb() {
		parents = new HashMap<String, Map<String, Integer>>();
		verbToRelations = new HashMap<String, Map<String, Double>>();
		typeToRelations = new HashMap<String, Map<String, Double>>();
		relationPriors = new HashMap<String, Double>();
		readHierarchy();
		readRelationPriors();
		readMapping();
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
		Map<Integer, Map<TokenSpan, Map<String, Double>>> billCats = new HashMap<Integer, Map<TokenSpan, Map<String, Double>>>();
        for (Triple<TokenSpan, String, Double> triple : document.getTokenSpanAnnotationConfidences(AnnotationTypeNLPCat.NELL_CATEGORY)) {
                        TokenSpan first = triple.getFirst();
                        int sentenceIndex = first.getSentenceIndex();
                        TokenSpan billCat = new TokenSpan(document, sentenceIndex, first.getStartTokenIndex(), first.getEndTokenIndex());
                        Map<TokenSpan, Map<String, Double>> tempMap = billCats.get(sentenceIndex);
                        if (tempMap == null) tempMap = new HashMap<TokenSpan, Map<String, Double>>();
                        Map<String, Double> cats = tempMap.get(billCat);
                        if (cats == null) cats = new HashMap<String, Double>();
                        cats.put(triple.getSecond(), triple.getThird());
                        tempMap.put(billCat, cats);
                        billCats.put(sentenceIndex, tempMap);
		}
        
        Map<Integer, Map<TokenSpan, Map<String, Double>>> oocCats = new HashMap<Integer, Map<TokenSpan, Map<String, Double>>>();
        for (Triple<TokenSpan, TypedNP, Double> triple : document.getTokenSpanAnnotationConfidences(ContextlessNPCategorizer.OUTOFCONTEXT_NP_CATEGORIES)) {
                        TokenSpan first = triple.getFirst();
                        int sentenceIndex = first.getSentenceIndex();
                        TokenSpan billCat = new TokenSpan(document, sentenceIndex, first.getStartTokenIndex(), first.getEndTokenIndex());
                        Map<TokenSpan, Map<String, Double>> tempMap = oocCats.get(sentenceIndex);
                        if (tempMap == null) tempMap = new HashMap<TokenSpan, Map<String, Double>>();
                        Map<String, Double> cats = tempMap.get(billCat);
                        if (cats == null) cats = new HashMap<String, Double>();
                        cats.put(triple.getSecond().getType(), triple.getThird());
                        tempMap.put(billCat, cats);
                        oocCats.put(sentenceIndex, tempMap);
		}        
        
        
        for (Map.Entry<Integer, Map<TokenSpan, Map<String, Double>>> e : oocCats.entrySet()) {
        	int sentenceIndex = e.getKey();
        	Map<TokenSpan, Map<String, Double>> currentMap = billCats.get(sentenceIndex);
        	if (currentMap == null) currentMap = new HashMap<TokenSpan, Map<String, Double>>();
        	Map<TokenSpan, Map<String, Double>> tempMap = e.getValue();
        	for (Map.Entry<TokenSpan, Map<String, Double>> e2 : tempMap.entrySet()) {
        		TokenSpan ts = e2.getKey();
        		if (currentMap.get(ts) == null) {
            		Map<String, Double> catMap = e2.getValue();
            		currentMap.put(ts, catMap);
        		}
        	}
        	if (currentMap.size() > 0) billCats.put(sentenceIndex, currentMap);
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
										if (!tempAL.contains(tokenVerbPPSpan)) tempAL.add(tokenVerbPPSpan);
										ArrayList<Integer> tempALI = foundVerbPPObjects.get(verbPP);
										if (tempALI == null) tempALI = new ArrayList<Integer>();
										if (!tempALI.contains(depIndex)) tempALI.add(depIndex);
										foundVerbPPs.put(verbPP, tempAL);
										foundVerbPPObjects.put(verbPP, tempALI);
										break;
									}
								}
							} else if (type.startsWith("nsubj") || type.equalsIgnoreCase("agent")) {
								if (!foundSubjects.contains(depIndex)) foundSubjects.add(depIndex);
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
								if (!foundDirectObjects.contains(depIndex)) foundDirectObjects.add(depIndex);
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
							String lemmaV = foundVerbs.get(i);
							if (isPassive) lemmaV = "(passive) " + lemmaV.trim();
							TokenSpan tokenSpanV = foundVerbTokenSpans.get(i);
							Map<String, Double> found = new HashMap<String, Double>();
							int obj = foundVerbObjects.get(i);
							String objectStr = document.getToken(sentIndex, obj).getStr();
							
							for (int subj : foundSubjects) {
								String subjectStr = document.getToken(sentIndex, subj).getStr();
								
								Map<TokenSpan, Map<String, Double>> tempMap = billCats.get(sentIndex);
								ArrayList<String> catSubj = null; ArrayList<String> catObj = null;
								ArrayList<String> topCatSubj = null; ArrayList<String> topCatObj = null;
								if (tempMap != null) {
									for (Map.Entry<TokenSpan, Map<String, Double>> e : tempMap.entrySet()) {
										TokenSpan ct = e.getKey();
										if (ct.containsToken(sentIndex, subj)) {
											catSubj = new ArrayList<String>();
											double maxVal = -1;
											for (Map.Entry<String, Double> ent : e.getValue().entrySet()) {
												if (!catSubj.contains(ent.getKey())) catSubj.add(ent.getKey());
												double val = ent.getValue();
												if (val > maxVal) {
													if (topCatSubj == null) topCatSubj = new ArrayList<String>();
													maxVal = val;
													topCatSubj.clear();
													if (!topCatSubj.contains(ent.getKey())) topCatSubj.add(ent.getKey());
												} else if (val == maxVal) {
													if (!topCatSubj.contains(ent.getKey())) topCatSubj.add(ent.getKey());
												}
											}
										}
										if (ct.containsToken(sentIndex, obj)) {
											catObj = new ArrayList<String>();
											double maxVal = -1;
											for (Map.Entry<String, Double> ent : e.getValue().entrySet()) {
												if (!catObj.contains(ent.getKey())) catObj.add(ent.getKey());
												double val = ent.getValue();
												if (val > maxVal) {
													if (topCatObj == null) topCatObj = new ArrayList<String>();
													maxVal = val;
													topCatObj.clear();
													if (!topCatObj.contains(ent.getKey())) topCatObj.add(ent.getKey());
												} else if (val == maxVal) {
													if (!topCatObj.contains(ent.getKey())) topCatObj.add(ent.getKey());
												}
											}
										}
									}
									
								}
								//System.out.println(catSubj.toString() + "\t" + catObj.toString() + "\t" + lemmaV);
								if (catSubj != null || catObj != null) {
									if (verbToRelations.get(lemmaV) != null) {
										Map<String, Double> candidateRelations = verbToRelations.get(lemmaV);
										if (catSubj != null && catObj != null) { // if both verb, subj and obj are found
											Map<String, Double> chosenCandidates = checkRelation(candidateRelations, catSubj, catObj);
											/*Map<String, Double> chosen = new HashMap<String, Double>();
											for (Map.Entry<String, Boolean> e : chosenCandidates.entrySet()) {
												chosen.put(e.getKey(), candidateRelations.get(e.getKey()));
											}
											Vector<String> sorted = sortMap(chosen);
											found.put(subjectStr + "\t" + catSubj.toString() +"||"+chosen.toString() + "\t" + sorted.toString()+"||"+objectStr + "\t" + catObj.toString(), 1.0);
											*/if (chosenCandidates.size() > 0) {
												Vector<String> sorted = sortMap(chosenCandidates);
												ArrayList<String> written = new ArrayList<String>();
												for (int j = sorted.size() - 1; j >=0; j--) {
													String relation = sorted.get(j);
													double conf = candidateRelations.get(relation);
													boolean domainFirst = relation.startsWith("dom::");
													String rrelation = relation.split(" ")[1].trim();
													if (!written.contains(rrelation)) {
														if (domainFirst) {
															found.put(subjectStr+"||"+rrelation+"||"+objectStr, conf);													
														} else {
															found.put(objectStr+"||"+rrelation+"||"+subjectStr, conf);
														}
														written.add(rrelation);
													}
												}
											}
											/*for (Map.Entry<String, Double> entry : candidateRelations.entrySet()) {
												String relation = entry.getKey().split(" ")[1];
												String chosen = checkRelation(relation, catSubj, catObj, true, true);
												double conf = entry.getValue();	
												
												if (chosen != null) found.put(subjectStr+"||"+chosen+"||"+objectStr, conf);
												String domain = entry.getKey().split(" ")[0];
												String range = entry.getKey().split(" ")[2];
												if (hasCategory(catSubj, domain) && hasCategory(catObj, range)) {
													//System.out.println(lemmaSubj + "(" + catSubj + ")" + "\t" + lemmaV + "\t" + lemmaObj + "(" + catObj + ")");
													verbs.add(new Triple<TokenSpan, String, Double>
														(new TokenSpan(document, sentIndex, tokenSpanV.getStartTokenIndex(), tokenSpanV.getEndTokenIndex()), 
																"svo:" +relation, conf));
													found.put("svo:"+chosen, conf);
													
												}
											}*/
										} else if (catSubj != null) { // if verb and subj are found
											/*for (Map.Entry<String, Double> entry : candidateRelations.entrySet()) {
												//String domain = entry.getKey().split(" ")[0];
												String relation = entry.getKey().split(" ")[1];
												double conf = entry.getValue();
												String chosen = checkRelation(relation, catSubj, catObj, true, false);
												if (chosen != null) found.put(subjectStr+"||"+chosen+"||"+objectStr, conf);
												if (hasCategory(catSubj, domain)) {
													verbs.add(new Triple<TokenSpan, String, Double>
														(new TokenSpan(document, sentIndex, tokenSpanV.getStartTokenIndex(), tokenSpanV.getEndTokenIndex()), 
																"sv:" + entry.getKey() + ":" + relation, conf));
													found.put("sv:" + relation, conf);
												}
											}*/
										} else if (catObj != null) { // if verb and obj are found
											/*for (Map.Entry<String, Double> entry : candidateRelations.entrySet()) {
												String relation = entry.getKey().split(" ")[1];
												//String range = entry.getKey().split(" ")[2];
												double conf = entry.getValue();
												String chosen = checkRelation(relation, catSubj, catObj, false, true);
												if (chosen != null) found.put(subjectStr+"||"+chosen+"||"+objectStr, conf);
												if (hasCategory(catObj, range)) {
													verbs.add(new Triple<TokenSpan, String, Double>
														(new TokenSpan(document, sentIndex, tokenSpanV.getStartTokenIndex(), tokenSpanV.getEndTokenIndex()), 
																"vo:" + relation, conf));
													found.put("vo:" + relation, conf);
												}
											}*/
										}
									} else if (topCatSubj != null && topCatObj != null) { // not found in verb-relation mapping, use the types to choose candidate relations
										/*for (String s : topCatSubj) {
											for (String o : topCatObj) {
												String type = s + ":::" + o;
												Map<String, Double> candidateRels = typeToRelations.get(type);
												if (candidateRels != null) {
													Vector<String> sortedCandidates = sortMap(candidateRels);
													ArrayList<Double> values = new ArrayList<Double>();
													for (int j = 0; j < sortedCandidates.size(); j++) {
														String relation = sortedCandidates.get(j);
														double conf = candidateRels.get(relation);
														if (!values.contains(conf)) values.add(conf);
														if (values.size() > 3) break;
														found.put(subjectStr+"||"+relation+"||"+objectStr, conf);
													}
													Iterator<Map.Entry<String, Double>> it1 = candidateRels.entrySet().iterator();
													while (it1.hasNext()) {
														Map.Entry<String, Double> ent1 = it1.next();
														String relation = ent1.getKey();
														double conf = ent1.getValue();
														found.put(subjectStr+"||"+relation+"||"+objectStr, conf);
													}
												} else {
													type = o + ":::" + s;
													candidateRels = typeToRelations.get(type);
													if (candidateRels != null) {
														Vector<String> sortedCandidates = sortMap(candidateRels);
														ArrayList<Double> values = new ArrayList<Double>();
														for (int j = 0; j < sortedCandidates.size(); j++) {
															String relation = sortedCandidates.get(j);
															double conf = candidateRels.get(relation);
															if (!values.contains(conf)) values.add(conf);
															if (values.size() > 3) break;
															found.put(objectStr+"||"+relation+"||"+subjectStr, conf);
														}
													}
												}
											}
										}*/
									}
								} 
							}
							if (found.size() == 0) { // if no subject/object; just use the verb without types
								/*if (verbToRelations.get(lemmaV) != null) {
									for (int subj : foundSubjects) {
										String subjectStr = document.getToken(sentIndex, subj).getStr();
										Map<String, Double> candidateRelations = verbToRelations.get(lemmaV);
										for (Map.Entry<String, Double> entry : candidateRelations.entrySet()) {
											String relation = entry.getKey().split(" ")[1];
											double conf = entry.getValue();
											verbs.add(new Triple<TokenSpan, String, Double>
											(new TokenSpan(document, sentIndex, tokenSpanV.getStartTokenIndex(), tokenSpanV.getEndTokenIndex()), 
													subjectStr+"||"+relation+"||"+objectStr, conf));
										}
									}
									
								}*/
							} else {
								for (Map.Entry<String, Double> e : found.entrySet()) {
									String relation = e.getKey();
									double conf = e.getValue();
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
		return verbs;
	}

	
	private Map<String, Double> checkRelation(
			Map<String, Double> candidateRelations, ArrayList<String> catSubj,
			ArrayList<String> catObj) {
		Map<String, Double> result = new HashMap<String, Double>();
		Map<String, Double> chosen = new HashMap<String, Double>();
		for (Map.Entry<String, Double> e : candidateRelations.entrySet()) {
			String relation = e.getKey();
			double conf = e.getValue();
			String domain = relation.split(" ")[0].trim().split("::")[1].trim();
			String range = relation.split(" ")[2].trim().split("::")[1].trim();
			boolean acceptedDomain = false; boolean acceptedRange = false;
			for (String cat : catSubj) {
				if (cat.equalsIgnoreCase(domain)) {
					acceptedDomain = true;
					break;
				}
				Map<String, Integer> parent = parents.get(cat);
				if (parent != null) {
					if (parent.get(domain) != null) {
						acceptedDomain = true;
						break;
					}
				}
			}
			for (String cat : catObj) {
				if (cat.equalsIgnoreCase(range)) {
					acceptedRange = true;
					break;
				}
				Map<String, Integer> parent = parents.get(cat);
				if (parent != null) {
					if (parent.get(range) != null) {
						acceptedRange = true;
						break;
					}
				}
			}
			if (acceptedDomain && acceptedRange) {
				chosen.put(relation, conf);
			}
		}
		
		/*Vector<String> sortedCandidates = sortMap(chosen);
		//ArrayList<Double> values = new ArrayList<Double>();
		for (int j = sortedCandidates.size() - 1; j >= 0; j--) {
			String relation = sortedCandidates.get(j);
			double conf = candidateRelations.get(relation);
			//if (!values.contains(conf)) values.add(conf);
			//if (values.size() > LIMIT) break;
			//boolean domainFirst = relation.startsWith("dom::");
			result.put(relation, conf);
		}*/
		return chosen;
	}

/*	private String checkRelation(String relation, ArrayList<String> catSubj,
			ArrayList<String> catObj, boolean subj, boolean obj) {
		double relVal = checkValue(relation, catSubj, catObj, subj, obj);
		if (relVal == 0) return relation;
		Vector<String> sorted = parentsSorted.get(relation);
		if (sorted != null) {
			for (String s : sorted) {
				double sVal = checkValue(s, catSubj, catObj, subj, obj);
				if (sVal >= 0 && sVal < relVal) return s;
			}
		}
		if (relVal >= 0) return relation;
		return null;
	}

	private double checkValue(String relation, ArrayList<String> catSubj,
			ArrayList<String> catObj, boolean subj, boolean obj) {
		double valDomain = -1;
		if (!subj) valDomain = 0;
		else {
			String domain = domains.get(relation);
			if (domain != null) {
				if (catSubj.contains(domain)) valDomain = 0;
				else {
					Vector<String> domainPs = parentsSorted.get(domain);
					if (domainPs != null) {
						for (String p : domainPs) {
							if (catSubj.contains(p)) {
								valDomain = (double) parents.get(domain).get(p);
								break;
							}
						}
					}
				}			
			}
		}

		double valRange = -1;
		if (!obj) valRange = 0;
		else {
			String range = ranges.get(relation);
			if (range != null) {
				if (catObj.contains(range)) valRange = 0;
				else {
					Vector<String> rangePs = parentsSorted.get(range);
					if (rangePs != null) {
						for (String p : rangePs) {
							if (catObj.contains(p)) {
								valRange = (double) parents.get(range).get(p);
								break;
							}
						}
					}
				}
			}
		}

		if (valDomain >= 0 && valRange >= 0) {
			double val = (valDomain + valRange) / (double) 2;
			return val;
		} 
		return -1;
	}*/

	private void readHierarchy() {
		Map<String, Map<String, Integer>> tempParents = new HashMap<String, Map<String, Integer>>();
		try {
			InputStream is = AnnotationVerb.class.getResourceAsStream("/NELL-hierarchy.txt");
			BufferedReader bfr = new BufferedReader(new InputStreamReader(is));
			String line, temp[];
			while ((line = bfr.readLine()) != null) {
				temp = line.split("\t");
				String child = temp[0].trim().replace("concept:", "");
				String parent = temp[2].trim().replace("concept:", "");
				
				Map<String, Integer> tempAL = tempParents.get(child);
				if (tempAL == null) tempAL = new HashMap<String, Integer>();
				tempAL.put(parent, 1);
				tempParents.put(child, tempAL);
			}
			bfr.close();
			
		} catch (IOException ioe) {
		    throw new RuntimeException(ioe);
		}
		
		for (Map.Entry<String, Map<String, Integer>> e : tempParents.entrySet()) {
			String node = e.getKey();
			Map<String, Integer> map = e.getValue();
			Map<String, Integer> processed = new HashMap<String, Integer>();
			ArrayList<String> toprocess = new ArrayList<String>();
			for (Map.Entry<String, Integer> e2 : map.entrySet()) {
				processed.put(e2.getKey(), e2.getValue());
				toprocess.add(e2.getKey());
			}
			while (toprocess.size() > 0) {
				String current = toprocess.remove(0);
				int val = processed.get(current);
				Map<String, Integer> map2 = tempParents.get(current);
				if (map2 != null) {
					for (Map.Entry<String, Integer> e3 : map2.entrySet()) {
						String candidate = e3.getKey();
						if (processed.get(candidate) == null) {
							processed.put(candidate, val+1);
							toprocess.add(candidate);
						}
					}
				}
			}
			parents.put(node, processed);
		}
	}

	private void readRelationPriors() {
		try {
			InputStream is = AnnotationVerb.class.getResourceAsStream("/priorProbRelation.txt");
			BufferedReader bfr = new BufferedReader(new InputStreamReader(is));
			String line, temp[];
			while ((line = bfr.readLine()) != null) {
				temp = line.split("\t");
				String relation = temp[1].trim().replace("concept:", "");
				double conf = Double.parseDouble(temp[2].trim());
				relationPriors.put(relation, conf);
			}
			bfr.close();
		} catch (IOException ioe) {
		    throw new RuntimeException(ioe);
		}
	}

	private void readMapping() {
		try {
			InputStream is = AnnotationVerb.class.getResourceAsStream("/naive-bayes-em-mapping-filtered-with-type-checking-and-prior-and-cps-with-decay-3-knee-verbs.txt");
			BufferedReader bfr = new BufferedReader(new InputStreamReader(is));
			String line;
			String temp[], relations[], verbs[], domain, range, relation, verb, type;
			double score;
			
			while ((line = bfr.readLine()) != null) {
				temp = line.split("\t");
				relations = temp[0].split(" ");
				domain = relations[0].trim().replace("concept:", "");
				range = relations[2].trim().replace("concept:", "");
				relation = relations[1].trim().replace("concept:", "");
				if (relationPriors.get(relation) != null) {
					type = domain + ":::" + range;
					Map<String, Double> relationMap = typeToRelations.get(type);
					if (relationMap == null) relationMap = new HashMap<String, Double>();
					relationMap.put(relation, relationPriors.get(relation));
					typeToRelations.put(type, relationMap);
				}
				temp = temp[1].split(";");
				for (String t : temp) {
					verbs = t.split(",");
					verb = verbs[0].trim();
					verb = verb.replace("(also in CPL)", "").trim();
					score = Double.parseDouble(verbs[1].trim());
					String rrelation;
					if (verb.startsWith("arg1")) {
						rrelation = "dom::" + domain + " " + relation + " " + "ran::" + range;
						verb = verb.replace("arg1 ", "").trim();
						verb = verb.replace(" arg2", "").trim();
					} else {
						rrelation = "ran::" + range + " " + relation + " " + "dom::" + domain;
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
	
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public static Vector<String> sortMap(Map<String, Double> map) {
		Map<String, Integer> sortedResults = sortByValue(map);
		List sortedQueries = new LinkedList(sortedResults.entrySet());
		Vector<String> sortedYears = new Vector<String>();
		
		for (Iterator itt = sortedQueries.iterator(); itt.hasNext();) {
			Map.Entry entry = (Map.Entry) itt.next();
			sortedYears.add((String) entry.getKey());
		}
		return sortedYears;
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	private static Map sortByValue(Map<String, Double> map) {
		List list = new LinkedList(map.entrySet());
		Collections.sort(list, new Comparator() {
			public int compare(Object o1, Object o2) {
				return ((Comparable) ((Map.Entry) (o1)).getValue()).compareTo(((Map.Entry) (o2)).getValue());
			}
		});
		Map result = new LinkedHashMap();
		for (Iterator it = list.iterator(); it.hasNext();) {
		     Map.Entry entry = (Map.Entry)it.next();
		     result.put(entry.getKey(), entry.getValue());
		     }
		return result;
	}
}