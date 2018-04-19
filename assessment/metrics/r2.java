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
public class r2 {
	public static double[] maioresValores;
	public static double[] menoresValores;
	
	public static void main(String[] args) throws IOException{
		int numExec=30;
		int objectiveNumber=-1;
		if(args.length == 0){
			System.err.println("uso: r2 front1 ... frontN");
			System.exit(1);
		}
// 		String problem=determineProblem(args);
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
//						System.out.println(calcular(lerReal(args[1]), dados));
						front[cont]=new Front();
						front[cont].solutions=new ArrayList<double[]>(dados);
						cont++;
					}
						
					dados.clear();
				}
				linha = br.readLine();
			}
			if(dados.size() > 0){//quando terminou o arquivo, processa os dados da ultima fronteira
				//System.out.println(calcular(lerReal(args[1]), dados));
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
		ArrayList<double[]> frontReal=new ArrayList<double[]>(lerReal("assessment/metrics/pareto/REF_"+objectiveNumber));
// // 		for(int s=0;s<frontReal.size();s++){//solucoes dentro do front real
// // 			double[] solTemp=frontReal.get(s);
// // 			for(int o=0;o<objectiveNumber;o++){ //objetivos da solucao
// // 				if(maioresValores[o]<solTemp[o])
// // 					maioresValores[o]=solTemp[o];
// // 				if(menoresValores[o]>solTemp[o])
// // 					menoresValores[o]=solTemp[o];
// // 			}
// // 		}
		
		//calcula tudo
		double[][] resultados = new double[file.size()][numExec];
// 		ArrayList<double[]> frontRealN = normalizeTudo(frontReal);
		for(int fe=0;fe<file.size();fe++){
			Front[] frontTemp = file.get(fe);
			for(int f=0;f<frontTemp.length;f++){ //fronts dentro do arquivo
				resultados[fe][f]=calcular(frontReal, normalizeTudo(frontTemp[f].solutions));
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
	
		static double max(double[] V, double[] A){
		double max=(Double.MAX_VALUE*-1);
		for(int j=0;j<V.length;j++){
			//double valor=V[j]*Math.abs(ponto[j]-A[j]);
			//double valor=V[j]*Math.abs(0-A[j]); //using the origin as reference point
			double valor1=V[j];
			double valor2=A[j];
			double valor = valor1*Math.abs(valor2);
			
			if( valor > max)
				max=valor;
		}
		//System.out.println("max: "+max);
		return max;
	}
	static double min(double[] V, ArrayList<double[]> fronteira){
		double min=Double.MAX_VALUE;
		for(int i=0;i<fronteira.size();i++){
			double valor=max(V, fronteira.get(i));
			if(valor < min)
				min=valor;
		}
		//System.out.println("min: "+min);
		return min;
	}
	
	public static double calcular(ArrayList<double[]> PFtrue, ArrayList<double[]> fronteira) {
		double R2=0;
		if(fronteira!=null){
			for(int i = 0; i<PFtrue.size(); i++){
				double[] pontoPFTrue = PFtrue.get(i);
				R2+=min(pontoPFTrue, fronteira);
				//System.out.println("R2: "+R2);
			}
			R2*=(1.0/PFtrue.size());
		} else{
			System.err.println("Erro no calculo do R2: Vetores de peso nao carregados.");
			System.exit(0);
			return 0;
		}	
		return R2;
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
	
// 	public static String determineProblem(String[] args){
// 		String problem="";
// 		String problemTest="";
// 		for(int i=0;i<args.length;i++){
// 			if(args[i].toLowerCase().contains("dtlz1"))
// 				problem="dtlz1";
// 			if(args[i].toLowerCase().contains("dtlz2"))
// 				problem="dtlz2";
// 			if(args[i].toLowerCase().contains("dtlz3"))
// 				problem="dtlz2";
// 			if(args[i].toLowerCase().contains("dtlz4"))
// 				problem="dtlz2";
// 			if(args[i].toLowerCase().contains("dtlz5"))
// 				problem="dtlz5";
// 			if(args[i].toLowerCase().contains("dtlz6"))
// 				problem="dtlz5";
// 			if(args[i].toLowerCase().contains("dtlz7"))
// 				problem="dtlz7";
// 			if(args[i].toLowerCase().contains("wfg1"))
// 				problem="wfg1";
// 			if(args[i].toLowerCase().contains("wfg2"))
// 				problem="wfg2";
// 			if(args[i].toLowerCase().contains("wfg3"))
// 				problem="wfg3";
// 			if(args[i].toLowerCase().contains("wfg4"))
// 				problem="wfg4";
// 			if(args[i].toLowerCase().contains("wfg5"))
// 				problem="wfg5";
// 			if(args[i].toLowerCase().contains("wfg6"))
// 				problem="wfg6";
// 			if(args[i].toLowerCase().contains("wfg7"))
// 				problem="wfg7";
// 			if(args[i].toLowerCase().contains("wfg8"))
// 				problem="wfg8";
// 			if(args[i].toLowerCase().contains("wfg9"))
// 				problem="wfg9";
// 			if(problemTest != "" && problemTest != problem){
// 				System.err.println("\nErro ao determinar o problema!!!\n");
// 				System.exit(1);
// 			}else
// 				problemTest=problem;
// 		}
// 		return problem;
// 	}
	
	public static double normalize(double valor, double min, double max){
		//return valor;
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
