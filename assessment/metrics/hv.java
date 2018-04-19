import java.lang.Math;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Arrays;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.FileWriter; 
import java.io.InputStreamReader;
import java.io.File;
import java.util.Random;

class Front{
	public ArrayList<double[]> solutions = new ArrayList<double[]>();
}
public class hv {
	public static double[] maioresValores;
	public static double[] menoresValores;
	public static double refPoint;
	
	public static void main(String[] args) throws IOException{
		int numExec=30;
		int objectiveNumber=-1;
		if(args.length == 0){
			System.err.println("uso: hv front1 ... frontN");
			System.exit(1);
		}
		String problem=determineProblem(args);
		String tempname=(problem+"-");
		ArrayList<Front[]> file = new ArrayList<Front[]>();
		
		for(int f=0;f<args.length;f++){
			//tempname+=args[f].split("/")[args[f].split("/").length-1];
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
		//ArrayList<double[]> frontReal=new ArrayList<double[]>(lerReal("pareto/DTLZ"+problem+"_"+objectiveNumber));
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
		//System.out.println(calcular(frontReal));
		//System.exit(1);
		//double hvReal=calcular(frontReal);
		for(int fe=0;fe<file.size();fe++){
			Front[] frontTemp = file.get(fe);
			for(int f=0;f<frontTemp.length;f++){ //fronts dentro do arquivo
				tempname=(tempname.split("-")[0])+"-"+(fe+"*"+f+"*"+objectiveNumber+"*"+System.nanoTime());
				resultados[fe][f]=calcular(frontTemp[f].solutions, tempname );
			}
		}
		//mostra os resultados
		for(int j=0;j<numExec;j++){
			for(int i=0;i<resultados.length;i++){
				//System.out.print((hvReal-resultados[i][j])+" ");
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
	
	public static double calcular(ArrayList<double[]> front, String tempname) throws IOException{
		double retorno=0;
		
		try{
			File file = new File("/tmp/temp"+tempname+".txt");
			FileWriter entrada = new FileWriter(file);
			entrada.write("#\n");
			for(int i=0;i<front.size();i++){
				for(int j=0;j<front.get(i).length;j++){
					entrada.write(normalize(front.get(i)[j], menoresValores[j], maioresValores[j])+" ");
				}
				entrada.write("\n");
			}
			entrada.write("#");
			entrada.flush();
			entrada.close();
			//Thread.sleep(1000);
			
			while(!file.exists()){
				Thread.sleep(1000);
			}
			
			if(!file.exists())
				System.err.println("File is missing before!");
	
			String ref="";
			refPoint=1.01;
			for(int i=0;i<front.get(0).length;i++)
				ref+=" 1.01";
				
			String comando="";
			if(front.get(0).length <= 10){
				comando+="assessment/metrics/hv/wfg /tmp/temp"+tempname+".txt"+ref+" | head -1";
				
			Process p =  Runtime.getRuntime().exec(new String[]{"/bin/sh", "-c", comando});
			p.waitFor();
			BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
			String line;
			while ( (line=br.readLine()) != null) {
				if(line.split("=").length == 2)
					retorno=Double.parseDouble(line.split("=")[1]);
				else{
					System.err.println("HV ERROR! ("+tempname+") - "+line);
					System.exit(1);
				}
			}
			
			if(!file.exists())
				System.err.println("File is missing after!");
			
			if(!file.delete())
				System.err.println("Temp file delete operation failed.");
			}else{
// 				comando+="java -Xmx2G -cp assessment/metrics/hv_mc Main /tmp/temp"+tempname+".txt"+ref+" | head -1";
				//comando+="assessment/metrics/hv_mc/a.out /tmp/temp"+tempname+".txt"+ref+" | head -1";
				retorno=calculaHipervolumeAmostragem(front);
			}
				
		}catch (Exception e){	e.printStackTrace();}
		return retorno;
	}
	
	
// 	public static double calcular(ArrayList<double[]> front) {
// 		double gdp = 0;
// 		if(fronteira!=null){
// 			
// 			double soma = 0;
// 			//Percorre todos os pontos do conjunto de aproximacao
// 			for(int i = 0; i<PFtrue.size(); i++){
// 				double[] pontoPFTrue = PFtrue.get(i);
// 				double menor_distancia = menorDistanciaEuclidiana(pontoPFTrue, fronteira); 
// 				soma+=  menor_distancia*menor_distancia;
// 				
// 			}
// 			
// 			//gd = Math.sqrt(soma)/(double) PFtrue.size(); //gd
// 			gdp = Math.sqrt( (1.0/PFtrue.size()) * soma ); //gdp
// 			
// 		} else{
// 			System.err.println("Erro no calculo do GD: Fronteira de Pareto nao carregada.");
// 			System.exit(0);
// 			return 0;
// 		}	
// 		return gdp;
// 	}

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


	//Calcula o hipervolume por Monte Carlo
	//Obs1. : Supõe que os dados estão normalizados seguindo o padrão do PISA, ou seja,
	//no intervalo [1.0,2.1]
	public static double calculaHipervolumeAmostragem(ArrayList<double[]> dados) throws IOException{
		Random random = new Random(System.nanoTime());
		int tamanhoAmostra=1000000;
	
		int dimensao = dados.get(0).length;
		int contaAcertos = 0;
		double gerado[] = new double[dimensao];
		// Modificar a linha abaixo caso os dados não estejam normalizados
		// ou estejam normalizados em um intervalo de tamanho diferente
		// 1.1 = 2.1 - 1.0, ou seja, limite superior - limite inferior
		
		//no meu caso 1.0 = 1 - 0
		
		// Se os dados tiverem diferentes limites para os diferentes objetivos
		// calcular o volumeTotal pelo produtório das diferenças dos limites....
		double volumeTotal = Math.pow(refPoint, dimensao);
		boolean dominado = true;


		for (int i = 0; i < tamanhoAmostra; i++) {
			for (int j = 0; j < dimensao; j++) {
				gerado[j] = refPoint * random.nextDouble();
			}

			dominado = false;

			//Verifica se o ponto é dominado ou não...
			for (int k = 0; ((k < dados.size()) && (!dominado)); k++) {
				double[] temp = dados.get(k);
				boolean dominadoTemp = true;

				for (int d = 0; d < dimensao; d++) {
					if (temp[d] > gerado[d]) {
						dominadoTemp = false;
					}
				}

				if (dominadoTemp) {
					dominado = true;
				}
			}

			if (dominado) {
				contaAcertos++;//conta o número de pontos dominados...
			}

		}

		//O hipervolume é proporcional a quantidade de pontos dominados e ao volume da área amostrada...
		return (double) contaAcertos / (double) tamanhoAmostra * volumeTotal;

	}
}