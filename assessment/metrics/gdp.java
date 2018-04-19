import java.lang.Math;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Arrays;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

class Front{
	public ArrayList<double[]> solutions = new ArrayList<double[]>();
}
public class gdp {
	public static double[] maioresValores;
	public static double[] menoresValores;
	
	public static void main(String[] args) throws IOException{
		int numExec=30;
		int objectiveNumber=-1;
		if(args.length == 0){
			System.err.println("uso: gdp front1 ... frontN");
			System.exit(1);
		}
		String problem=determineProblem(args);
		ArrayList<Front[]> file = new ArrayList<Front[]>();
		
		for(int f=0;f<args.length;f++){
			FileReader fr = new FileReader(args[f]);
			BufferedReader br = new BufferedReader(fr);
			ArrayList<double[]> dados = new ArrayList<double[]>();
			Front[] front = new Front[numExec];
			int cont = 0;

			String linha = br.readLine();
			while (linha != null) { //nao chegou no fim do arquivo
				String[] objs = linha.split("\\s+"); //separa por espacos
				if(objs.length==0) //se nao deu
					objs = linha.split("\t");//tenta separar por tab
				if ((objs.length != 1 || !(objs[0].equals(""))) && !objs[0].equals("#") ) { //se tiver um caractere, nao vazio e diferente de #
					objectiveNumber=objs.length;
					double[] d = new double[objectiveNumber];
					for (int i = 0; i < objectiveNumber; i++) {
						d[i] = Double.parseDouble(objs[i].replace(",",".")); // converte a linha para numeros tanto com ponto como com virgula de separador
					}
					dados.add(d);
				} else { //se chegar na linha em branco ou invalida
					if(dados.size() > 0){
						//System.out.println(calcular(dados, lerReal(args[1])));
						front[cont]=new Front();
						front[cont].solutions=new ArrayList<double[]>(dados);
						cont++;
					}
						
					dados.clear();
				}
				linha = br.readLine();
			}
			if(dados.size() > 0){//quando terminou o arquivo, processa os dados da ultima fronteira
				//System.out.println(calcular(dados, lerReal(args[1])));
				front[cont]=new Front();
				front[cont].solutions=new ArrayList<double[]>(dados);
				cont++;
			}
			file.add(front);
		}
		
		//System.out.println("arquivos: "+file.size());
		//inicializa os maiores e menores valores
		maioresValores = new double[objectiveNumber];
		menoresValores = new double[objectiveNumber];
		for(int o=0;o<objectiveNumber;o++){ //objetivos da solucao
			maioresValores[o]=Double.MAX_VALUE*-1;
			menoresValores[o]=Double.MAX_VALUE;
		}
		//atualiza os maiores e menores valores
		for(int fe=0;fe<file.size();fe++){
			Front[] frontTemp = file.get(fe);
			for(int f=0;f<frontTemp.length;f++){ //fronts dentro do arquivo
				for(int s=0;s<frontTemp[f].solutions.size();s++){//solucoes dentro do front
					double[] solTemp=frontTemp[f].solutions.get(s);
					for(int o=0;o<objectiveNumber;o++){ //objetivos da solucao
						if(maioresValores[o]<solTemp[o])
							maioresValores[o]=solTemp[o];
						if(menoresValores[o]>solTemp[o])
							menoresValores[o]=solTemp[o];
					}
				}
			}
		}
		//le o front real e atualiza os maiores e menores valores (pra incluir 0 como origem na maioria)
		ArrayList<double[]> frontReal=new ArrayList<double[]>(lerReal("assessment/metrics/pareto/"+problem.toUpperCase()+"_"+objectiveNumber));
		for(int s=0;s<frontReal.size();s++){//solucoes dentro do front real
			double[] solTemp=frontReal.get(s);
			for(int o=0;o<objectiveNumber;o++){ //objetivos da solucao
				if(maioresValores[o]<solTemp[o])
					maioresValores[o]=solTemp[o];
				if(menoresValores[o]>solTemp[o])
					menoresValores[o]=solTemp[o];
			}
		}
		
		//calcula tudo
		double[][] resultados = new double[file.size()][numExec];
		ArrayList<double[]> frontRealN = normalizeTudo(frontReal);
		for(int fe=0;fe<file.size();fe++){
			Front[] frontTemp = file.get(fe);
			for(int f=0;f<frontTemp.length;f++){ //fronts dentro do arquivo
				resultados[fe][f]=calcular(normalizeTudo(frontTemp[f].solutions), frontRealN );
			}
		}
		//mostra os resultados
		for(int j=0;j<numExec;j++){
			for(int i=0;i<resultados.length;i++){
				System.out.print(resultados[i][j]+" ");
			}
			System.out.println();
		}
		
// 		System.out.println("maiores: ");
// 		for(int o=0;o<objectiveNumber;o++){ //objetivos da solucao
// 			System.out.print(maioresValores[o]+" ");
// 		}
// 		System.out.println("\nmenores: ");
// 		for(int o=0;o<objectiveNumber;o++){ //objetivos da solucao
// 			System.out.print(menoresValores[o]+" ");
// 		}
// 		System.out.println();
		
	}
	
	public static double calcular(ArrayList<double[]> PFtrue, ArrayList<double[]> fronteira) {
		double gdp = 0;
		if(fronteira!=null){
			
			double soma = 0;
			//Percorre todos os pontos do conjunto de aproximacao
			for(int i = 0; i<PFtrue.size(); i++){
				double[] pontoPFTrue = PFtrue.get(i);
				double menor_distancia = menorDistanciaEuclidiana(pontoPFTrue, fronteira); 
				soma+=  menor_distancia*menor_distancia;
				
			}
			
			//gd = Math.sqrt(soma)/(double) PFtrue.size(); //gd
			gdp = Math.sqrt( (1.0/PFtrue.size()) * soma ); //gdp
			
		} else{
			System.err.println("Erro no calculo do GD: Fronteira de Pareto nao carregada.");
			System.exit(0);
			return 0;
		}	
		return gdp;
	}


	protected static double menorDistanciaEuclidiana(double[] ponto, ArrayList<double[]> fronteira) {
		double menor_distancia = Double.MAX_VALUE;
		//Obtem a menor dist�ncia entre um ponto do conjunto de aproxima��o e um ponto da fronteira de pareto real
		for(int j = 0; j<fronteira.size();j++){
			double[] ponto2 = fronteira.get(j);
			menor_distancia = Math.min(distanciaEuclidiana(ponto, ponto2), menor_distancia);
		}
		return menor_distancia;
	}
	public static ArrayList<double[]> lerReal(String arquivo) throws IOException{
		FileReader fr = new FileReader(arquivo);
		BufferedReader br = new BufferedReader(fr);
		ArrayList<double[]> dados = new ArrayList<double[]>();
		int cont = 0;

		String linha = br.readLine();
		while (linha != null) {
			String[] objs = linha.split("\\s+");
			if(objs.length==0)
				objs = linha.split("\t");
			if ((objs.length != 1 || !(objs[0].equals(""))) && !objs[0].equals("#") ) {
				double[] d = new double[objs.length];
				for (int i = 0; i < objs.length; i++) {
					d[i] = Double.parseDouble(objs[i].replace(",","."));
				}
				dados.add(d);
			} else {
				if(dados.size() > 0)
					//System.out.println(calcular(dados));
					return dados;

				dados.clear();
			}
			linha = br.readLine();
		}
//		if(dados.size() > 0)
			return dados;
	}
	public static double distanciaEuclidiana(double[] vetor1, double[] vetor2){
		double soma = 0;
		for (int i = 0; i < vetor1.length; i++) {
			double valor1=vetor1[i];
			double valor2=vetor2[i];
			//soma += Math.pow(vetor1[i]-vetor2[i],2);
			soma += (valor1-valor2)*(valor1-valor2);
		}
		return Math.sqrt(soma);
	}
	
	public static String determineProblem(String[] args){
		String problem="";
		String problemTest="";
		for(int i=0;i<args.length;i++){
			if(args[i].toLowerCase().contains("dtlz1"))
				problem="dtlz1";
			if(args[i].toLowerCase().contains("dtlz2"))
				problem="dtlz2";
			if(args[i].toLowerCase().contains("dtlz3"))
				problem="dtlz2";
			if(args[i].toLowerCase().contains("dtlz4"))
				problem="dtlz2";
			if(args[i].toLowerCase().contains("dtlz5"))
				problem="dtlz5";
			if(args[i].toLowerCase().contains("dtlz6"))
				problem="dtlz5";
			if(args[i].toLowerCase().contains("dtlz7"))
				problem="dtlz7";
			if(args[i].toLowerCase().contains("wfg1"))
				problem="wfg1";
			if(args[i].toLowerCase().contains("wfg2"))
				problem="wfg2";
			if(args[i].toLowerCase().contains("wfg3"))
				problem="wfg3";
			if(args[i].toLowerCase().contains("wfg4"))
				problem="wfg4";
			if(args[i].toLowerCase().contains("wfg5"))
				problem="wfg5";
			if(args[i].toLowerCase().contains("wfg6"))
				problem="wfg6";
			if(args[i].toLowerCase().contains("wfg7"))
				problem="wfg7";
			if(args[i].toLowerCase().contains("wfg8"))
				problem="wfg8";
			if(args[i].toLowerCase().contains("wfg9"))
				problem="wfg9";
			if(args[i].toLowerCase().contains("uf1"))
				problem="uf1";
			if(args[i].toLowerCase().contains("uf2"))
				problem="uf2";
			if(args[i].toLowerCase().contains("uf3"))
				problem="uf3";
			if(args[i].toLowerCase().contains("uf4"))
				problem="uf4";
			if(args[i].toLowerCase().contains("uf5"))
				problem="uf5";
			if(args[i].toLowerCase().contains("uf6"))
				problem="uf6";
			if(args[i].toLowerCase().contains("uf7"))
				problem="uf7";
			if(args[i].toLowerCase().contains("uf8"))
				problem="uf8";
			if(args[i].toLowerCase().contains("uf9"))
				problem="uf9";
			if(args[i].toLowerCase().contains("uf10"))
				problem="uf10";
			if(problemTest != "" && problemTest != problem){
				System.err.println("\nErro ao determinar o problema!!!\n");
				System.exit(1);
			}else
				problemTest=problem;
		}
		return problem;
	}
	
	public static double normalize(double valor, double min, double max){
		if(max-min > 0){
			return ((valor-min)/(max-min));
		}else{
			if(valor > 0)
				System.err.println("min: "+min+", max: "+max+", valor: "+valor+"\n");
			return (valor-min);//((valor-min)/(max-min));
		}
	}
	public static ArrayList<double[]> normalizeTudo(ArrayList<double[]> entrada){
		ArrayList<double[]> saida = new ArrayList<double[]>();
	
		for(int s=0;s<entrada.size();s++){//solucoes dentro do front de entrada
			double[] sd = new double[entrada.get(s).length];
			for(int o=0;o<sd.length;o++){ //objetivos da solucao
				sd[o]=normalize(entrada.get(s)[o], menoresValores[o], maioresValores[o]);
			}
			saida.add(sd);
		}
		return saida;
	}
}
