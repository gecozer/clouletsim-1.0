/*
 * Title:        CloudSim Toolkit
 * Description:  CloudSim (Cloud Simulation) Toolkit for Modeling and Simulation of Clouds
 * Licence:      GPL - http://www.gnu.org/copyleft/gpl.html
 *
 * Copyright (c) 2009-2012, The University of Melbourne, Australia
 */

package org.cloudbus.cloudsim;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.cloudbus.cloudsim.core.CloudSim;
import org.cloudbus.cloudsim.core.CloudSimTags;
import org.cloudbus.cloudsim.core.SimEntity;
import org.cloudbus.cloudsim.core.SimEvent;
import org.cloudbus.cloudsim.lists.CloudletList;
import org.cloudbus.cloudsim.lists.VmList;

/**
 * DatacentreBroker represents a broker acting on behalf of a user. It hides VM management, as vm
 * creation, sumbission of cloudlets to this VMs and destruction of VMs.
 * 
 * @author Rodrigo N. Calheiros
 * @author Anton Beloglazov
 * @since CloudSim Toolkit 1.0
 */
public class DatacenterBroker extends SimEntity {

	/** The vm list. */
	protected List<? extends Vm> vmList;

	/** The vms created list. */
	protected List<? extends Vm> vmsCreatedList;

	/** The cloudlet list. */
	protected List<? extends Cloudlet> cloudletList;

	/** The cloudlet submitted list. */
	protected List<? extends Cloudlet> cloudletSubmittedList;

	/** The cloudlet received list. */
	protected List<? extends Cloudlet> cloudletReceivedList;

	/** The cloudlets submitted. */
	protected int cloudletsSubmitted;

	/** The vms requested. */
	protected int vmsRequested;

	/** The vms acks. */
	protected int vmsAcks;

	/** The vms destroyed. */
	protected int vmsDestroyed;

	/** The datacenter ids list. */
	protected List<Integer> datacenterIdsList;

	/** The datacenter requested ids list. */
	protected List<Integer> datacenterRequestedIdsList;

	/** The vms to datacenters map. */
	protected Map<Integer, Integer> vmsToDatacentersMap;

	/** The datacenter characteristics list. */
	protected Map<Integer, DatacenterCharacteristics> datacenterCharacteristicsList;

	/**
	 * Created a new DatacenterBroker object.
	 * 
	 * @param name name to be associated with this entity (as required by Sim_entity class from
	 *            simjava package)
	 * @throws Exception the exception
	 * @pre name != null
	 * @post $none
	 */
	public DatacenterBroker(String name) throws Exception {
		super(name);

		setVmList(new ArrayList<Vm>());
		setVmsCreatedList(new ArrayList<Vm>());
		setCloudletList(new ArrayList<Cloudlet>());
		setCloudletSubmittedList(new ArrayList<Cloudlet>());
		setCloudletReceivedList(new ArrayList<Cloudlet>());

		cloudletsSubmitted = 0;
		setVmsRequested(0);
		setVmsAcks(0);
		setVmsDestroyed(0);

		setDatacenterIdsList(new LinkedList<Integer>());
		setDatacenterRequestedIdsList(new ArrayList<Integer>());
		setVmsToDatacentersMap(new HashMap<Integer, Integer>());
		setDatacenterCharacteristicsList(new HashMap<Integer, DatacenterCharacteristics>());
	}

	/**
	 * This method is used to send to the broker the list with virtual machines that must be
	 * created.
	 * 
	 * @param list the list
	 * @pre list !=null
	 * @post $none
	 */
	public void submitVmList(List<? extends Vm> list) {
		getVmList().addAll(list);
	}

	/**
	 * This method is used to send to the broker the list of cloudlets.
	 * 
	 * @param list the list
	 * @pre list !=null
	 * @post $none
	 */
	public void submitCloudletList(List<? extends Cloudlet> list) {
		getCloudletList().addAll(list);
	}

	/**
	 * Specifies that a given cloudlet must run in a specific virtual machine.
	 * 
	 * @param cloudletId ID of the cloudlet being bount to a vm
	 * @param vmId the vm id
	 * @pre cloudletId > 0
	 * @pre id > 0
	 * @post $none
	 */
	public void bindCloudletToVm(int cloudletId, int vmId) {
		CloudletList.getById(getCloudletList(), cloudletId).setVmId(vmId);
	}
	//˳���㷨
	public void bindCloudletsToVmsSimple(){
		int vmNum = vmList.size();
		int cloudletNum = cloudletList.size();
		int idx = 0;
		for(int i=0;i<cloudletNum;i++){
		cloudletList.get(i).setVmId(vmList.get(idx).getId());
		idx = (idx+1)%vmNum;
		}
	}
	//̰���㷨
	public void bindCloudletsToVmsTimeAwared(){
		int cloudletNum=cloudletList.size();
		int vmNum=vmList.size();
		//time[i][j] ��ʾ����i�������j�ϵ�ִ��ʱ��
		double[][] time=new double[cloudletNum][vmNum];
		//cloudletList��MI��������, vm��MIPS��������
		Collections.sort(cloudletList,new CloudletComparator());
		Collections.sort(vmList,new VmComparator());
		//////////For test//////////////////////////////////
		System.out.println("///////////For test///////////////");
		for(int i=0;i<cloudletNum;i++){
			System.out.print(cloudletList.get(i).getCloudletId()+":"+cloudletList.get
			(i).getCloudletLength()+" ");
			}
			System.out.println();
			for(int i=0;i<vmNum;i++){
			System.out.print(vmList.get(i).getId()+":"+vmList.get(i).getMips()+" ");
			}
			System.out.println();
			System.out.println("//////////////////////////////////");
			//////////////////////////////////////////////////////////////////////
			for(int i=0;i<cloudletNum;i++){
			for(int j=0;j<vmNum;j++){
			time[i][j]=
			(double)cloudletList.get(i).getCloudletLength()/vmList.get(j).getMips();
			System.out.print("time["+i+"]["+j+"]="+time[i][j]+" ");
			//For test
			}
			System.out.println(); //For test
			}
			double[] vmLoad=new double[vmNum];//��ĳ����������������ִ��ʱ��
			int[] vmTasks=new int[vmNum]; //��ĳ��Vm�����е���������
			double minLoad=0;//��¼��ǰ������䷽ʽ������ֵ
			int idx=0;//��¼��ǰ�������ŷ��䷽ʽ��Ӧ��������к�
			//��һ��cloudlet���������vm
			vmLoad[vmNum-1]=time[0][vmNum-1];
			vmTasks[vmNum-1]=1;
			cloudletList.get(0).setVmId(vmList.get(vmNum-1).getId());
			for(int i=1;i<cloudletNum;i++){
				minLoad=vmLoad[vmNum-1]+time[i][vmNum-1];
				idx=vmNum-1;
				for(int j=vmNum-2;j>=0;j--){
				//�����ǰ�����δ����������Ƚ��굱ǰ����������������Ƿ�����
					if(vmLoad[j]==0){
						if(minLoad>=time[i][j])idx=j;
							break;
					}
					if(minLoad>vmLoad[j]+time[i][j]){
						minLoad=vmLoad[j]+time[i][j];
						idx=j;
					}
					//�򵥵ĸ��ؾ���
					else if(minLoad==vmLoad[j]+time[i][j]&&vmTasks[j]<vmTasks[idx])
						idx=j;
				}
				vmLoad[idx]+=time[i][idx];
				vmTasks[idx]++;
				cloudletList.get(i).setVmId(vmList.get(idx).getId());
				System.out.print(i+"th "+"vmLoad["+idx+"]="+vmLoad[idx]+"minLoad="+minLoad);
				System.out.println();
			}
		}
	public void bindCloudletsToVmsACO(){
		//��ÿֻ����Я��n�����񣬽������������У���ʼʱ��n���������
		//���䵽������ϣ�ÿ������Ҫ���ľ���������һ����������ҵ�һ��������ͷ���һ������
		//ͬʱ������Ϣ�صĸ��£�ÿ���һ��һ����Ҫ����İ������¼������棬һ����Ϣ�صı�ʾ
		//������ֹ����������
		//Ϊ��ʹ����򵥣�ʹ�����ʱ����һ��һ������Ϊ·���û��ĺ�����׼
		//�����Ƽ��������������Դ��������������
		int cloudletNum=cloudletList.size();//cloudletList�����б�
		
		int n=vmList.size();//vmList�������Դ�б�
		//time[i][j] ��ʾ����i�������j�ϵ�ִ��ʱ��
		double[][] time=new double[cloudletNum][n];

		double Alpha=0.9;//������Ϣ����Ҫ�̶ȵĲ���
		double Beta=2;//��������ʽ������Ҫ�̶ȵĲ���
		double Rho=0.1;//��Ϣ������ϵ��
		double Rho_min=0.05;
		int NC_max=100;//���һ������
		int Q=50;//��Ϣ������ǿ��ϵ��
		int m=31;//���ϸ���
		int R_best[][]=new int[200][cloudletNum];//�������·��
		double L_best[]=new double[NC_max];//�������ʱ��
		double Eta[][]=new double[n][cloudletNum];;//EtaΪ��������
		double Tau[][]=new double[n][n];//TauΪ��Ϣ�ؾ���
		for(int i=0;i<NC_max;i++) L_best[i]=10000;
		
		int Tabu[][]=new int[m][cloudletNum];//�洢����¼������б�����ɣ����ɱ�
		for(int i=0;i<m;i++)
			for(int j=0;j<cloudletNum;j++)
				Tabu[i][j]=0;
		int it=1;
		double D[][]=new double[n][cloudletNum];//�������Ϊ�������ӵ���Դ
		List<Integer> s=new ArrayList<Integer>();
		List<Integer> s1=new ArrayList<Integer>();
		double min=10000;
		int position=0;
		
		double totalLength=0,totalInputfileSize=0;
		double q1=0.3,q2=0.3,q3=0.4;//��QosԼ����������ȣ���Ϊʱ�䣬�ɱ����Լ��ɿ����������棬���յ�Ȩ�ط���Ӧ���ɿͻ�����������
		int ii=0,jj=0;
		double tmp_min=1000,tmp_max=0;
		for(int i=0;i<vmList.size();i++){
			if(tmp_max<vmList.get(i).getMips()) {
				tmp_max=vmList.get(i).getMips();
				ii=i;
			}
			if(tmp_min>vmList.get(i).getMips()){
				tmp_min=vmList.get(i).getMips();
				jj=i;
			}
		}
		double completeTime_globalmin=0,completeTime_globalmax=0;
		double completeCost_globalmin=0,completeCost_globalmax=0;
		
		
		for(int k=0;k<cloudletList.size();k++){
			totalLength+=cloudletList.get(k).getCloudletLength();
			totalInputfileSize+=cloudletList.get(k).getCloudletFileSize();
		}
		completeTime_globalmin=totalLength/(cloudletList.size()*tmp_max)+totalInputfileSize/(cloudletList.size()*vmList.get(ii).getBw());
		completeTime_globalmax=totalLength/(cloudletList.size()*tmp_min)+totalInputfileSize/(cloudletList.size()*vmList.get(jj).getBw());
		completeCost_globalmin=completeTime_globalmin*1;
		completeCost_globalmax=completeTime_globalmin*1;
		
		for(int i2=0;i2<n;i2++)
			for(int j2=0;j2<cloudletNum;j2++){
				D[i2][j2]=cloudletList.get(j2).getCloudletLength()/vmList.get(i2).getMips()
					+cloudletList.get(j2).getCloudletFileSize()/vmList.get(i2).getBw();
				//Ӧ��������������ɵģ�һ���������С�����еľ��������ļ���С
				Eta[i2][j2]=1/D[i2][j2];
			}
		for(int i=0;i<n;i++)
			s1.add(i);
		while(it<=NC_max){
			//������������������������
			s.clear();
			for(int i=0;i<n;i++)
				for(int j=0;j<n;j++)
					Tau[i][j]=1;//��Ϣ�س�ʼ��Ϊ1
			for(int i=0;i<(m/n+1);i++){
				Collections.shuffle(s1);
				s.addAll(s1);
			}
			//System.out.println(s);
			int isVolocate=1;
			for(int i=0;i<m;i++) Tabu[i][0]=s.get(i);//��������ֲ���������ϣ�ÿֻ���ϵ�һ�����ڵ�λ��
			for(int i=1;i<cloudletNum;i++){//��tsp���ⲻһ�����ǣ�ÿֻ����ѭ���Ĵ���Ӧ�������������Ŀ�������ģ�����������ĸ���
				//���յ����ʱ��Ӧ����ÿ������������������ɵ����ʱ��
				for(int j=0;j<m;j++){
					//������¼�Ѿ������������ĸ����������������Ѿ�������ʱ�ô�����
					List<Integer> visited=new ArrayList<Integer>();
					for(int x=0;x<i;x++) visited.add(Tabu[j][x]);//��¼�ѷ���������
					List<Integer> J=new ArrayList<Integer>();//�������ɴ����������
					List<Double> P=new ArrayList<Double>();//�������ɴ�����������ĸ���
					boolean temp=false;
					
					for(int i2=0;i2<n;i2++){
						int tag=0;
						if(i>n-1)//���������Ѿ���������
							//����Ӧ�������е�����������Ա�����
						{
							if(!temp)
								J.clear();
							J.add(i2);
							temp=true;
							
						}
						else{
							for(int k=0;k<visited.size();k++){
								if(visited.get(k)==i2) tag=1;//δ���������
								
							}
							if(tag==0) J.add(i2);
						}
						
					}
						//����ÿ�������������ĸ���
					for(int k=0;k<J.size();k++)
						//����ط�Ӧ�ó����������ӣ�eta��beta
						P.add(Math.pow(Tau[visited.get(visited.size()-1)][J.get(k)],Alpha)
								*Math.pow(Eta[visited.get(visited.size()-1)][J.get(k)],Beta));
					double sum=0;
					int tovm=0;
					List<Double> tempP=new ArrayList<Double>();
					for(int k=0;k<P.size();k++) sum+=P.get(k);
					for(int k=0;k<P.size();k++) P.set(k, P.get(k)/sum);
					tempP.addAll(P);
					//���Ų������̶ĵķ�ʽѡ����һ�����
					for(int i2=0;i2<P.size();i2++)
					{
						sum=0;
						for(int j2=0;j2<=i2;j2++)
						{
							sum+=tempP.get(j2);
						}
						P.set(i2, sum);
					}
					for(int i2=0;i2<P.size();i2++){
						if(P.get(i2)>Math.random()){
							tovm=J.get(i2);//��¼��һѡ��������������������Ǽ�¼��J�е�
							break;
						}
					}
					Tabu[j][i]=tovm;
					
				}

				
				
				//����ط����Խ��оֲ���Ϣ�صĸ���,ͬʱ���϶�ά��Լ��
				double completeTime_localmin=0,completeTime_localmax=0;
				double completeCost_localmin=0,completeCost_localmax=0;
				
				
				for(int i2=0;i2<m;i2++){
					double completeTime=0,completeCost=0,completeRes=0;
					double resTime=0,resCost=0,resRes=0;
					double Qos=0;//Լ������ֵ
					double time2[]=new double[n];
					totalLength=0;
					totalInputfileSize=0;
					
					if(i<5){
						for(int k=0;k<=i;k++){
						if(completeTime<cloudletList.get(k).getCloudletLength()/vmList.get(Tabu[i2][k]).getMips()
								+cloudletList.get(k).getCloudletFileSize()/vmList.get(Tabu[i2][k]).getBw())
						completeTime=cloudletList.get(k).getCloudletLength()/vmList.get(Tabu[i2][k]).getMips()
							+cloudletList.get(k).getCloudletFileSize()/vmList.get(Tabu[i2][k]).getBw();
						}
					}
					else{
						for(int k2=0;k2<n;k2++)
							for(int k=0;k<=i;k++){
								if(Tabu[i2][k]==k2) time2[k2]+=cloudletList.get(k).getCloudletLength()/vmList.get(k2).getMips()
								 	+cloudletList.get(k).getCloudletFileSize()/vmList.get(Tabu[i2][k]).getBw();
							}
						for(int k2=0;k2<n;k2++)
							if(completeTime<time2[k2]) completeTime=time2[k2];
					}
					for(int k=0;k<=i;k++){
							totalLength+=cloudletList.get(k).getCloudletLength();
							totalInputfileSize+=cloudletList.get(k).getCloudletFileSize();
							
							completeCost+=cloudletList.get(k).getCloudletLength()/vmList.get(Tabu[i2][k]).getMips()*1
								+cloudletList.get(k).getCloudletFileSize()/vmList.get(Tabu[i2][k]).getBw()*1;
					    
					    	completeRes+=(cloudletList.get(k).getCloudletLength()/vmList.get(Tabu[i2][k]).getMips())/500;
						}
					completeTime_localmin=totalLength/((i+1)*tmp_max)+totalInputfileSize/((i+1)*vmList.get(0).getBw());
					completeTime_localmax=totalLength/((i+1)*tmp_min)+totalInputfileSize/((i+1)*vmList.get(0).getBw());
					completeCost_localmin=completeTime_localmin*1;
					completeCost_localmax=completeTime_localmax*1;
					
					
					resTime=(completeTime-completeTime_localmin)/(completeTime_localmax-completeTime_localmin);
					resCost=(completeCost-completeCost_localmin)/(completeCost_localmax-completeCost_localmin);
					resRes=completeRes/m;
					Qos=q1*resTime+q2*resCost+q3*resRes;
					for(int j2=0;j2<=i-1;j2++)
						Tau[Tabu[i2][j2]][Tabu[i2][j2+1]]=(1-Rho)*Tau[Tabu[i2][j2]][Tabu[i2][j2+1]]+Q/Qos;
				}
				
				
				
				
				isVolocate+=1;
			}
			
			if(it>2){
				for(int i2=0;i2<cloudletNum;i2++)
					Tabu[1][i2]=R_best[it-1][i2];
			}
			//�ҳ������������
			int L[]=new int[m];
			//ÿ����������������õ���ʱ��,ע��Ӧ���ǲ��е�
			
			for(int i2=0;i2<m;i2++){
				List<Integer> R=new ArrayList<Integer>();
				int tmpL[]=new int[n];
				int tmp=0;
				for(int j2=0;j2<cloudletNum;j2++) R.add(Tabu[i2][j2]);
				for(int k=0;k<n;k++){//����ÿ�������ִ������������ʱ��
					for(int p=0;p<R.size();p++)
						if(R.get(p)==k){
							tmpL[k]+=D[k][p];
						}
				}
				for(int k=0;k<n;k++)
					if(tmp<tmpL[k]) tmp=tmpL[k];
				L[i2]=tmp;
			}
			
			min=10000;
			position=0;
			for(int i2=0;i2<m;i2++){
				
				
				if(min>L[i2]) {
					min=L[i2];
					position=i2;//�������
				}
			}
			L_best[it-1]=min;//ÿ�������ʱ��
			//����ط�����������Ӧ�Ļӷ�ϵ�����޸ģ�ͨ���Ƚ�ÿ��֮��Ĳ���
			//��Ϣ�ظߵ�Ӧ���ǻӷ����ٶ��Ǹ����
			if(it-1>=4){
				if(L_best[it-1]>L_best[it-5]){
					if(Rho>=Rho_min)
						Rho=0.95*Rho;
					else
						Rho=Rho_min;
				}
			}
			
			for(int i2=0;i2<cloudletNum;i2++) R_best[it-1][i2]=Tabu[position][i2];
			it=it+1;
			//������Ϣ��
			/*int update_tau[][]=new int[n][n];
			for(int i2=0;i2<m;i2++)
				for(int j2=0;j2<n-1;j2++)
				{
					update_tau[Tabu[i2][j2]][Tabu[i2][j2+1]]=update_tau[Tabu[i2][j2]][Tabu[i2][j2+1]]+Q/L[i2];
				}
			for(int i2=0;i2<n;i2++)
				for(int j2=0;j2<n;j2++)
				{
					Tau[i2][j2]=(1-Rho)*update_tau[i2][j2]+update_tau[i2][j2];//��Ϣ���лӷ�
				}*/
			//��ֲ�����Ϣ�ظ��������ж�άԼ����ȫ����Ϣ�ظ���
			
			for(int i2=0;i2<m;i2++){
				double completeTime=0,completeCost=0,completeRes=0;
				double resTime=0,resCost=0,resRes=0;
				double Qos=0;//Լ������ֵ
				double time2[]=new double[n];
				for(int k2=0;k2<n;k2++)
					for(int k=0;k<cloudletList.size();k++){
						if(Tabu[i2][k]==k2) time2[k2]+=cloudletList.get(k).getCloudletLength()/vmList.get(k2).getMips()
								 	+cloudletList.get(k).getCloudletFileSize()/vmList.get(Tabu[i2][k]).getBw();
						
						for(int k3=0;k3<n;k3++)
							if(completeTime<time2[k3]) completeTime=time2[k3];
					}
				for(int k=0;k<cloudletList.size();k++){
					completeCost+=cloudletList.get(k).getCloudletLength()/vmList.get(Tabu[i2][k]).getMips()*1
						+cloudletList.get(k).getCloudletFileSize()/vmList.get(Tabu[i2][k]).getBw()*1;
				    completeRes+=(cloudletList.get(k).getCloudletLength()/vmList.get(Tabu[i2][k]).getMips())/500;
				}
				completeCost=completeTime*1;
				resTime=(completeTime-completeTime_globalmin)/(completeTime_globalmax-completeTime_globalmin);
				resCost=(completeCost-completeCost_globalmin)/(completeCost_globalmax-completeCost_globalmin);
				resRes=completeRes/m;
				Qos=q1*resTime+q2*resCost+q3*resRes;
				for(int j2=0;j2<cloudletList.size()-1;j2++)
					Tau[Tabu[i2][j2]][Tabu[i2][j2+1]]=(1-Rho)*Tau[Tabu[i2][j2]][Tabu[i2][j2+1]]+Q/Qos;
			}
		
		
			
			for(int i2=0;i2<m;i2++)
				for(int j2=0;j2<cloudletNum;j2++)
				{
					Tabu[i2][j2]=0;//���ɱ�����
				}
			Rho=Math.random();
		}
		//��ӡ�������ѵ���������䷽ʽ���Լ�ʱ��
		min=10000;
		position=0;
		for(int i2=0;i2<NC_max;i2++){
			if(min>L_best[i2]) {
				min=L_best[i2];
				position=i2;//��Ѵ���
			}
		}
		System.out.println("���ʱ�䣺"+min);
		System.out.print("��������У�");
		for(int i2=0;i2<cloudletNum;i2++){
			System.out.print(R_best[position][i2]+"  ");
			cloudletList.get(i2).setVmId(vmList.get(R_best[position][i2]).getId());
		}
		
		System.out.println();
		double cost=0;
		for(int i=0;i<cloudletList.size();i++){
			cost+=cloudletList.get(i).getCloudletLength()/vmList.get(R_best[position][i]).getMips()*0.9
				+cloudletList.get(i).getCloudletFileSize()/vmList.get(R_best[position][i]).getBw()*0.5;
		}
		System.out.println("�ɱ�:"+cost);
		int load[]=new int[n];
		double sysload=0,avgload=cloudletList.size()/n;
		int tmp=0;
		for(int i=0;i<n;i++){
			tmp=0;
			for(int j=0;j<cloudletList.size();j++)
				if(R_best[position][j]==i) tmp++;
			load[i]=tmp;
		}
		for(int i=0;i<n;i++){
			sysload+=Math.pow(load[i]-avgload, 2);
		}
		sysload=Math.pow(sysload, 0.5)/n;
		System.out.println("ϵͳ��Դ����:"+sysload);
		System.out.println("--------------------------------------");
	}
	public void bindCloudletsToVmsACONormal(){
		//��ÿֻ����Я��n�����񣬽������������У���ʼʱ��n���������
		//���䵽������ϣ�ÿ������Ҫ���ľ���������һ����������ҵ�һ��������ͷ���һ������
		//ͬʱ������Ϣ�صĸ��£�ÿ���һ��һ����Ҫ����İ������¼������棬һ����Ϣ�صı�ʾ
		//������ֹ����������
		//Ϊ��ʹ����򵥣�ʹ�����ʱ����һ��һ������Ϊ·���û��ĺ�����׼
		//�����Ƽ��������������Դ��������������
		int cloudletNum=cloudletList.size();//cloudletList�����б�
		
		int n=vmList.size();//vmList�������Դ�б�
		//time[i][j] ��ʾ����i�������j�ϵ�ִ��ʱ��
		double[][] time=new double[cloudletNum][n];

		double Alpha=0.9;//������Ϣ����Ҫ�̶ȵĲ���
		int Beta=2;//��������ʽ������Ҫ�̶ȵĲ���
		double Rho=0.1;//��Ϣ������ϵ��
		int NC_max=200;//���һ������
		int Q=100;//��Ϣ������ǿ��ϵ��
		int m=31;//���ϸ���
		int R_best[][]=new int[200][cloudletNum];//�������·��
		double L_best[]=new double[NC_max];//�������ʱ��
		double Eta[][]=new double[n][cloudletNum];;//EtaΪ��������
		double Tau[][]=new double[n][n];//TauΪ��Ϣ�ؾ���
		for(int i=0;i<NC_max;i++) L_best[i]=10000;
		for(int i=0;i<n;i++)
			for(int j=0;j<n;j++)
				Tau[i][j]=1;//��Ϣ�س�ʼ��Ϊ1
		int Tabu[][]=new int[m][cloudletNum];//�洢����¼������б�����ɣ����ɱ�
		for(int i=0;i<m;i++)
			for(int j=0;j<cloudletNum;j++)
				Tabu[i][j]=0;
		int it=1;
		double D[][]=new double[n][cloudletNum];//�������Ϊ�������ӵ���Դ
		List<Integer> s=new ArrayList<Integer>();
		List<Integer> s1=new ArrayList<Integer>();
		double min=10000;
		int position=0;
		for(int i2=0;i2<n;i2++)
			for(int j2=0;j2<cloudletNum;j2++){
				D[i2][j2]=cloudletList.get(j2).getCloudletLength()/vmList.get(i2).getMips()
					+cloudletList.get(j2).getCloudletFileSize()/vmList.get(i2).getBw();
				Eta[i2][j2]=1/D[i2][j2];
			}
		for(int i=0;i<n;i++)
			s1.add(i);
		while(it<=NC_max){
			//������������������������
			s.clear();
			for(int i=0;i<(m/n+1);i++){
				Collections.shuffle(s1);
				s.addAll(s1);
			}
			//System.out.println(s);
			for(int i=0;i<m;i++) Tabu[i][0]=s.get(i);//��������ֲ���������ϣ�ÿֻ���ϵ�һ�����ڵ�λ��
			for(int i=1;i<cloudletNum;i++){//��tsp���ⲻһ�����ǣ�ÿֻ����ѭ���Ĵ���Ӧ�������������Ŀ�������ģ�����������ĸ���
				//���յ����ʱ��Ӧ����ÿ������������������ɵ����ʱ��
				for(int j=0;j<m;j++){
					int isVolocate=1;//������¼�Ѿ������������ĸ����������������Ѿ�������ʱ�ô�����
					List<Integer> visited=new ArrayList<Integer>();
					for(int x=0;x<i;x++) visited.add(Tabu[j][x]);//��¼�ѷ���������
					List<Integer> J=new ArrayList<Integer>();//�������ɴ����������
					List<Double> P=new ArrayList<Double>();//�������ɴ�����������ĸ���
					boolean temp=false;
					
					for(int i2=0;i2<n;i2++){
						int tag=0;
						if(i>n-1)//���������Ѿ���������
							//����Ӧ�������е�����������Ա�����
						{
							if(!temp)
								J.clear();
							J.add(i2);
							temp=true;
							
						}
						else{
							for(int k=0;k<visited.size();k++){
								if(visited.get(k)==i2) tag=1;//δ���������
								
							}
							if(tag==0) J.add(i2);
						}
						
					}
						//����ÿ�������������ĸ���
					for(int k=0;k<J.size();k++)
						//����ط�Ӧ�ó����������ӣ�eta��beta
						P.add(Math.pow(Tau[visited.get(visited.size()-1)][J.get(k)],Alpha)
								*Math.pow(Eta[visited.get(visited.size()-1)][J.get(k)],Beta));
					double sum=0;
					int tovm=0;
					List<Double> tempP=new ArrayList<Double>();
					for(int k=0;k<P.size();k++) sum+=P.get(k);
					for(int k=0;k<P.size();k++) P.set(k, P.get(k)/sum);
					tempP.addAll(P);
					//���Ų������̶ĵķ�ʽѡ����һ�����
					for(int i2=0;i2<P.size();i2++)
					{
						sum=0;
						for(int j2=0;j2<=i2;j2++)
						{
							sum+=tempP.get(j2);
						}
						P.set(i2, sum);
					}
					for(int i2=0;i2<P.size();i2++){
						if(P.get(i2)>Math.random()){
							tovm=J.get(i2);//��¼��һѡ��������������������Ǽ�¼��J�е�
							break;
						}
					}
					Tabu[j][i]=tovm;
					isVolocate+=1;
				}
			}
			if(it>2){
				for(int i2=0;i2<cloudletNum;i2++)
					Tabu[1][i2]=R_best[it-1][i2];
			}
			//�ҳ������������
			int L[]=new int[m];
			//ÿ����������������õ���ʱ��,ע��Ӧ���ǲ��е�
			
			for(int i2=0;i2<m;i2++){
				List<Integer> R=new ArrayList<Integer>();
				int tmpL[]=new int[n];
				int tmp=0;
				for(int j2=0;j2<cloudletNum;j2++) R.add(Tabu[i2][j2]);
				for(int k=0;k<n;k++){//����ÿ�������ִ������������ʱ��
					for(int p=0;p<R.size();p++)
						if(R.get(p)==k){
							tmpL[k]+=D[k][p];
						}
				}
				for(int k=0;k<n;k++)
					if(tmp<tmpL[k]) tmp=tmpL[k];
				L[i2]=tmp;
			}
			
			min=10000;
			position=0;
			for(int i2=0;i2<m;i2++){
				
				
				if(min>L[i2]) {
					min=L[i2];
					position=i2;//�������
				}
			}
			L_best[it-1]=min;//ÿ�������ʱ��
			for(int i2=0;i2<cloudletNum;i2++) R_best[it-1][i2]=Tabu[position][i2];
			it=it+1;
			//������Ϣ��
			int update_tau[][]=new int[n][n];
			for(int i2=0;i2<m;i2++)
				for(int j2=0;j2<n-1;j2++)
				{
					update_tau[Tabu[i2][j2]][Tabu[i2][j2+1]]=update_tau[Tabu[i2][j2]][Tabu[i2][j2+1]]+Q/L[i2];
				}
			for(int i2=0;i2<n;i2++)
				for(int j2=0;j2<n;j2++)
				{
					Tau[i2][j2]=(1-Rho)*update_tau[i2][j2]+update_tau[i2][j2];//��Ϣ���лӷ�
				}
			
			for(int i2=0;i2<m;i2++)
				for(int j2=0;j2<cloudletNum;j2++)
				{
					Tabu[i2][j2]=0;//���ɱ�����
				}
		}
		//��ӡ�������ѵ���������䷽ʽ���Լ�ʱ��
		min=10000;
		position=0;
		for(int i2=0;i2<NC_max;i2++){
			if(min>L_best[i2]) {
				min=L_best[i2];
				position=i2;//��Ѵ���
			}
		}
		System.out.println("���ʱ�䣺"+min);
		System.out.print("��������У�");
		for(int i2=0;i2<cloudletNum;i2++){
			System.out.print(R_best[position][i2]+"  ");
			cloudletList.get(i2).setVmId(vmList.get(R_best[position][i2]).getId());
		}
		
		System.out.println();
		double cost=0;
		for(int i=0;i<cloudletList.size();i++){
			cost+=cloudletList.get(i).getCloudletLength()/vmList.get(R_best[position][i]).getMips()*0.9
				+cloudletList.get(i).getCloudletFileSize()/vmList.get(R_best[position][i]).getBw()*0.5;
		}
		System.out.println("�ɱ�:"+cost);
		int load[]=new int[n];
		double sysload=0,avgload=cloudletList.size()/n;
		int tmp=0;
		for(int i=0;i<n;i++){
			tmp=0;
			for(int j=0;j<cloudletList.size();j++)
				if(R_best[position][j]==i) tmp++;
			load[i]=tmp;
		}
		for(int i=0;i<n;i++){
			sysload+=Math.pow(load[i]-avgload, 2);
		}
		sysload=Math.pow(sysload, 0.5)/n;
		System.out.println("ϵͳ��Դ����:"+sysload);
		System.out.println("--------------------------------------");
	}
	public void bindCloudletsToVmsSA(){//ģ���˻��㷨
		int cloudletNum=cloudletList.size();//cloudletList�����б�
		int n=vmList.size();
		double D[][]=new double[n][cloudletNum];//����������Ӧ����������ķѵ�ʱ��
		for(int i2=0;i2<n;i2++)
			for(int j2=0;j2<cloudletNum;j2++){
				D[i2][j2]=cloudletList.get(j2).getCloudletLength()/vmList.get(i2).getMips()
					+cloudletList.get(j2).getCloudletFileSize()/vmList.get(i2).getBw();
		}
		List<Integer> s1=new ArrayList<Integer>();
		List<Integer> s=new ArrayList<Integer>();
		List<Integer> s0=new ArrayList<Integer>();
		double Sum=10000;//�Ϻõ�����ʱ��
		double tmp=0;
		for(int i=0;i<n;i++)
			s1.add(i);
		for(int j=0;j<10;j++){//ͨ�����ַ������һ���Ϻõĳ�ʼ��
			
			s.clear();
			for(int i=0;i<cloudletNum/n;i++){
				Collections.shuffle(s1);
				
				s.addAll(s1);
			}
			Collections.shuffle(s1);
			for(int i=0;i<cloudletNum-(cloudletNum/n)*n;i++)
				s.add(s1.get(i));//��ó�ʼ���������
			//����ʱ�䣬ע���ǲ��е�
			
			tmp=0;
			double tmpL[]=new double[n];
			for(int k=0;k<n;k++){//����ÿ�������ִ������������ʱ��
				for(int p=0;p<s.size();p++)
					if(s.get(p)==k){
						tmpL[k]+=D[k][p];
					}
			}
			for(int k=0;k<n;k++)
				if(tmp<tmpL[k]) tmp=tmpL[k];
			if(tmp<Sum) {
				Sum=tmp;
				s0.clear();
				s0.addAll(s);
			}
		}
		//����ģ���˻�����ϵ��
		double e=Math.pow(0.1, 30);
		int L=200;
		double at=0.999;
		double T=1;
		for(int i=0;i<L;i++){	//�˻����
			int min=0,max=0;
			Random rand=new Random();
			double df=0;//���ۺ���ֵ
			double before=0,after=0;
			double tmpL2[]=new double[n];
			int tmp2=0;
			int p=rand.nextInt(cloudletNum-1),q=rand.nextInt(cloudletNum-1);
			List<Integer> s0_tmp=new ArrayList<Integer>();
			s0_tmp.addAll(s0);
			do{
				if(p==q) q=rand.nextInt(cloudletNum-1);
			}while(p==q);
			if(p<q){ max=q;min=p;}
			else{max=p;min=q;};
			for(int k=0;k<n;k++){//����ÿ�������ִ������������ʱ��
				for(int p1=0;p1<s0.size();p1++)
					if(s0.get(p1)==k){
						tmpL2[k]+=D[k][p1];
					}
			}
			for(int k=0;k<n;k++)
				if(before<tmpL2[k]) before=tmpL2[k];
			//����
			for(int k=0;k<max-min+1;k++){
				tmp2=s0.get(min+k);
				s0_tmp.set(min+k, s0.get(max-k));
				s0_tmp.set(max-k, tmp2);
			}
			
			for(int k=0;k<n;k++){//����ÿ�������ִ������������ʱ��
				tmpL2[k]=0;
				for(int p1=0;p1<s0_tmp.size();p1++)
					if(s0_tmp.get(p1)==k){
						tmpL2[k]+=D[k][p1];
					}
			}
			for(int k=0;k<n;k++)
				if(after<tmpL2[k]) after=tmpL2[k];
			df=after-before;
			if(df<0){ 
				Sum=Sum+df;
				s0.clear();
				s0.addAll(s0_tmp);
			}
			else{
				if(Math.pow(Math.E, -df/T)>Math.random())  { Sum=Sum+df;s0.clear();s0.addAll(s0_tmp); };
			}
			T=T*at;
			if(T<e) break;
			s0_tmp.clear();
		}
		double tmpL3[]=new double[n];
		tmp=0;
		for(int i2=0;i2<cloudletNum;i2++){
			System.out.print(s0.get(i2)+"-");
		}
		System.out.println();
		double cost=0;
		for(int i=0;i<cloudletList.size();i++){
			cost+=cloudletList.get(i).getCloudletLength()/vmList.get(s0.get(i)).getMips()*0.9
				+cloudletList.get(i).getCloudletFileSize()/vmList.get(s0.get(i)).getBw()*0.5;
		}
		System.out.println("�ɱ�:"+cost);
		int load[]=new int[n];
		double sysload=0,avgload=cloudletList.size()/n;
		int tmpp=0;
		for(int i=0;i<n;i++){
			tmpp=0;
			for(int j=0;j<cloudletList.size();j++)
				if(s0.get(j)==i) tmpp++;
			load[i]=tmpp;
		}
		for(int i=0;i<n;i++){
			sysload+=Math.pow(load[i]-avgload, 2);
		}
		sysload=Math.pow(sysload, 0.5)/n;
		System.out.println("ϵͳ��Դ����:"+sysload);
		for(int k=0;k<n;k++){//����ÿ�������ִ������������ʱ��
			for(int p=0;p<s0.size();p++)
				if(s0.get(p)==k){
					tmpL3[k]+=D[k][p];
				}
		}
		for(int k=0;k<n;k++)
			if(tmp<tmpL3[k]) tmp=tmpL3[k];
		System.out.println("�����ѵ�ʱ�䣺"+tmp);
		for(int i2=0;i2<cloudletNum;i2++){
			cloudletList.get(i2).setVmId(s0.get(i2));
		}
		
		System.out.println();
		System.out.println("--------------------------------------");
	}
			//Cloudlet����MI��������
	private class CloudletComparator implements Comparator<Cloudlet>{
			public int compare(Cloudlet cl1,Cloudlet cl2){
			return (int)(cl2.getCloudletLength()-cl1.getCloudletLength());
			}
			}
			//Vm����MIPS��������
	private class VmComparator implements Comparator<Vm>{
			public int compare(Vm vm1,Vm vm2){
			return (int)(vm1.getMips()-vm2.getMips());
			}
			}
			/////////////////////////////////////////////////////////////////////////
	/**
	 * Processes events available for this Broker.
	 * 
	 * @param ev a SimEvent object
	 * @pre ev != null
	 * @post $none
	 */
	@Override
	public void processEvent(SimEvent ev) {
		switch (ev.getTag()) {
		// Resource characteristics request
			case CloudSimTags.RESOURCE_CHARACTERISTICS_REQUEST:
				processResourceCharacteristicsRequest(ev);
				break;
			// Resource characteristics answer
			case CloudSimTags.RESOURCE_CHARACTERISTICS:
				processResourceCharacteristics(ev);
				break;
			// VM Creation answer
			case CloudSimTags.VM_CREATE_ACK:
				processVmCreate(ev);
				break;
			// A finished cloudlet returned
			case CloudSimTags.CLOUDLET_RETURN:
				processCloudletReturn(ev);
				break;
			// if the simulation finishes
			case CloudSimTags.END_OF_SIMULATION:
				shutdownEntity();
				break;
			// other unknown tags are processed by this method
			default:
				processOtherEvent(ev);
				break;
		}
	}

	/**
	 * Process the return of a request for the characteristics of a PowerDatacenter.
	 * 
	 * @param ev a SimEvent object
	 * @pre ev != $null
	 * @post $none
	 */
	protected void processResourceCharacteristics(SimEvent ev) {
		DatacenterCharacteristics characteristics = (DatacenterCharacteristics) ev.getData();
		getDatacenterCharacteristicsList().put(characteristics.getId(), characteristics);

		if (getDatacenterCharacteristicsList().size() == getDatacenterIdsList().size()) {
			setDatacenterRequestedIdsList(new ArrayList<Integer>());
			createVmsInDatacenter(getDatacenterIdsList().get(0));
		}
	}

	/**
	 * Process a request for the characteristics of a PowerDatacenter.
	 * 
	 * @param ev a SimEvent object
	 * @pre ev != $null
	 * @post $none
	 */
	protected void processResourceCharacteristicsRequest(SimEvent ev) {
		setDatacenterIdsList(CloudSim.getCloudResourceList());
		setDatacenterCharacteristicsList(new HashMap<Integer, DatacenterCharacteristics>());

		Log.printLine(CloudSim.clock() + ": " + getName() + ": Cloud Resource List received with "
				+ getDatacenterIdsList().size() + " resource(s)");

		for (Integer datacenterId : getDatacenterIdsList()) {
			sendNow(datacenterId, CloudSimTags.RESOURCE_CHARACTERISTICS, getId());
		}
	}

	/**
	 * Process the ack received due to a request for VM creation.
	 * 
	 * @param ev a SimEvent object
	 * @pre ev != null
	 * @post $none
	 */
	protected void processVmCreate(SimEvent ev) {
		int[] data = (int[]) ev.getData();
		int datacenterId = data[0];
		int vmId = data[1];
		int result = data[2];

		if (result == CloudSimTags.TRUE) {
			getVmsToDatacentersMap().put(vmId, datacenterId);
			getVmsCreatedList().add(VmList.getById(getVmList(), vmId));
			Log.printLine(CloudSim.clock() + ": " + getName() + ": VM #" + vmId
					+ " has been created in Datacenter #" + datacenterId + ", Host #"
					+ VmList.getById(getVmsCreatedList(), vmId).getHost().getId());
		} else {
			Log.printLine(CloudSim.clock() + ": " + getName() + ": Creation of VM #" + vmId
					+ " failed in Datacenter #" + datacenterId);
		}

		incrementVmsAcks();

		// all the requested VMs have been created
		if (getVmsCreatedList().size() == getVmList().size() - getVmsDestroyed()) {
			submitCloudlets();
		} else {
			// all the acks received, but some VMs were not created
			if (getVmsRequested() == getVmsAcks()) {
				// find id of the next datacenter that has not been tried
				for (int nextDatacenterId : getDatacenterIdsList()) {
					if (!getDatacenterRequestedIdsList().contains(nextDatacenterId)) {
						createVmsInDatacenter(nextDatacenterId);
						return;
					}
				}

				// all datacenters already queried
				if (getVmsCreatedList().size() > 0) { // if some vm were created
					submitCloudlets();
				} else { // no vms created. abort
					Log.printLine(CloudSim.clock() + ": " + getName()
							+ ": none of the required VMs could be created. Aborting");
					finishExecution();
				}
			}
		}
	}

	/**
	 * Process a cloudlet return event.
	 * 
	 * @param ev a SimEvent object
	 * @pre ev != $null
	 * @post $none
	 */
	protected void processCloudletReturn(SimEvent ev) {
		Cloudlet cloudlet = (Cloudlet) ev.getData();
		getCloudletReceivedList().add(cloudlet);
		Log.printLine(CloudSim.clock() + ": " + getName() + ": Cloudlet " + cloudlet.getCloudletId()
				+ " received");
		cloudletsSubmitted--;
		if (getCloudletList().size() == 0 && cloudletsSubmitted == 0) { // all cloudlets executed
			Log.printLine(CloudSim.clock() + ": " + getName() + ": All Cloudlets executed. Finishing...");
			clearDatacenters();
			finishExecution();
		} else { // some cloudlets haven't finished yet
			if (getCloudletList().size() > 0 && cloudletsSubmitted == 0) {
				// all the cloudlets sent finished. It means that some bount
				// cloudlet is waiting its VM be created
				clearDatacenters();
				createVmsInDatacenter(0);
			}

		}
	}

	/**
	 * Overrides this method when making a new and different type of Broker. This method is called
	 * by {@link #body()} for incoming unknown tags.
	 * 
	 * @param ev a SimEvent object
	 * @pre ev != null
	 * @post $none
	 */
	protected void processOtherEvent(SimEvent ev) {
		if (ev == null) {
			Log.printLine(getName() + ".processOtherEvent(): " + "Error - an event is null.");
			return;
		}

		Log.printLine(getName() + ".processOtherEvent(): "
				+ "Error - event unknown by this DatacenterBroker.");
	}

	/**
	 * Create the virtual machines in a datacenter.
	 * 
	 * @param datacenterId Id of the chosen PowerDatacenter
	 * @pre $none
	 * @post $none
	 */
	protected void createVmsInDatacenter(int datacenterId) {
		// send as much vms as possible for this datacenter before trying the next one
		int requestedVms = 0;
		String datacenterName = CloudSim.getEntityName(datacenterId);
		for (Vm vm : getVmList()) {
			if (!getVmsToDatacentersMap().containsKey(vm.getId())) {
				Log.printLine(CloudSim.clock() + ": " + getName() + ": Trying to Create VM #" + vm.getId()
						+ " in " + datacenterName);
				sendNow(datacenterId, CloudSimTags.VM_CREATE_ACK, vm);
				requestedVms++;
			}
		}

		getDatacenterRequestedIdsList().add(datacenterId);

		setVmsRequested(requestedVms);
		setVmsAcks(0);
	}

	/**
	 * Submit cloudlets to the created VMs.
	 * 
	 * @pre $none
	 * @post $none
	 */
	protected void submitCloudlets() {
		int vmIndex = 0;
		for (Cloudlet cloudlet : getCloudletList()) {
			Vm vm;
			// if user didn't bind this cloudlet and it has not been executed yet
			if (cloudlet.getVmId() == -1) {
				vm = getVmsCreatedList().get(vmIndex);
			} else { // submit to the specific vm
				vm = VmList.getById(getVmsCreatedList(), cloudlet.getVmId());
				if (vm == null) { // vm was not created
					Log.printLine(CloudSim.clock() + ": " + getName() + ": Postponing execution of cloudlet "
							+ cloudlet.getCloudletId() + ": bount VM not available");
					continue;
				}
			}

			Log.printLine(CloudSim.clock() + ": " + getName() + ": Sending cloudlet "
					+ cloudlet.getCloudletId() + " to VM #" + vm.getId());
			cloudlet.setVmId(vm.getId());
			sendNow(getVmsToDatacentersMap().get(vm.getId()), CloudSimTags.CLOUDLET_SUBMIT, cloudlet);
			cloudletsSubmitted++;
			vmIndex = (vmIndex + 1) % getVmsCreatedList().size();
			getCloudletSubmittedList().add(cloudlet);
		}

		// remove submitted cloudlets from waiting list
		for (Cloudlet cloudlet : getCloudletSubmittedList()) {
			getCloudletList().remove(cloudlet);
		}
	}

	/**
	 * Destroy the virtual machines running in datacenters.
	 * 
	 * @pre $none
	 * @post $none
	 */
	protected void clearDatacenters() {
		for (Vm vm : getVmsCreatedList()) {
			Log.printLine(CloudSim.clock() + ": " + getName() + ": Destroying VM #" + vm.getId());
			sendNow(getVmsToDatacentersMap().get(vm.getId()), CloudSimTags.VM_DESTROY, vm);
		}

		getVmsCreatedList().clear();
	}

	/**
	 * Send an internal event communicating the end of the simulation.
	 * 
	 * @pre $none
	 * @post $none
	 */
	protected void finishExecution() {
		sendNow(getId(), CloudSimTags.END_OF_SIMULATION);
	}

	/*
	 * (non-Javadoc)
	 * @see cloudsim.core.SimEntity#shutdownEntity()
	 */
	@Override
	public void shutdownEntity() {
		Log.printLine(getName() + " is shutting down...");
	}

	/*
	 * (non-Javadoc)
	 * @see cloudsim.core.SimEntity#startEntity()
	 */
	@Override
	public void startEntity() {
		Log.printLine(getName() + " is starting...");
		schedule(getId(), 0, CloudSimTags.RESOURCE_CHARACTERISTICS_REQUEST);
	}

	/**
	 * Gets the vm list.
	 * 
	 * @param <T> the generic type
	 * @return the vm list
	 */
	@SuppressWarnings("unchecked")
	public <T extends Vm> List<T> getVmList() {
		return (List<T>) vmList;
	}

	/**
	 * Sets the vm list.
	 * 
	 * @param <T> the generic type
	 * @param vmList the new vm list
	 */
	protected <T extends Vm> void setVmList(List<T> vmList) {
		this.vmList = vmList;
	}

	/**
	 * Gets the cloudlet list.
	 * 
	 * @param <T> the generic type
	 * @return the cloudlet list
	 */
	@SuppressWarnings("unchecked")
	public <T extends Cloudlet> List<T> getCloudletList() {
		return (List<T>) cloudletList;
	}

	/**
	 * Sets the cloudlet list.
	 * 
	 * @param <T> the generic type
	 * @param cloudletList the new cloudlet list
	 */
	protected <T extends Cloudlet> void setCloudletList(List<T> cloudletList) {
		this.cloudletList = cloudletList;
	}

	/**
	 * Gets the cloudlet submitted list.
	 * 
	 * @param <T> the generic type
	 * @return the cloudlet submitted list
	 */
	@SuppressWarnings("unchecked")
	public <T extends Cloudlet> List<T> getCloudletSubmittedList() {
		return (List<T>) cloudletSubmittedList;
	}

	/**
	 * Sets the cloudlet submitted list.
	 * 
	 * @param <T> the generic type
	 * @param cloudletSubmittedList the new cloudlet submitted list
	 */
	protected <T extends Cloudlet> void setCloudletSubmittedList(List<T> cloudletSubmittedList) {
		this.cloudletSubmittedList = cloudletSubmittedList;
	}

	/**
	 * Gets the cloudlet received list.
	 * 
	 * @param <T> the generic type
	 * @return the cloudlet received list
	 */
	@SuppressWarnings("unchecked")
	public <T extends Cloudlet> List<T> getCloudletReceivedList() {
		return (List<T>) cloudletReceivedList;
	}

	/**
	 * Sets the cloudlet received list.
	 * 
	 * @param <T> the generic type
	 * @param cloudletReceivedList the new cloudlet received list
	 */
	protected <T extends Cloudlet> void setCloudletReceivedList(List<T> cloudletReceivedList) {
		this.cloudletReceivedList = cloudletReceivedList;
	}

	/**
	 * Gets the vm list.
	 * 
	 * @param <T> the generic type
	 * @return the vm list
	 */
	@SuppressWarnings("unchecked")
	public <T extends Vm> List<T> getVmsCreatedList() {
		return (List<T>) vmsCreatedList;
	}

	/**
	 * Sets the vm list.
	 * 
	 * @param <T> the generic type
	 * @param vmsCreatedList the vms created list
	 */
	protected <T extends Vm> void setVmsCreatedList(List<T> vmsCreatedList) {
		this.vmsCreatedList = vmsCreatedList;
	}

	/**
	 * Gets the vms requested.
	 * 
	 * @return the vms requested
	 */
	protected int getVmsRequested() {
		return vmsRequested;
	}

	/**
	 * Sets the vms requested.
	 * 
	 * @param vmsRequested the new vms requested
	 */
	protected void setVmsRequested(int vmsRequested) {
		this.vmsRequested = vmsRequested;
	}

	/**
	 * Gets the vms acks.
	 * 
	 * @return the vms acks
	 */
	protected int getVmsAcks() {
		return vmsAcks;
	}

	/**
	 * Sets the vms acks.
	 * 
	 * @param vmsAcks the new vms acks
	 */
	protected void setVmsAcks(int vmsAcks) {
		this.vmsAcks = vmsAcks;
	}

	/**
	 * Increment vms acks.
	 */
	protected void incrementVmsAcks() {
		vmsAcks++;
	}

	/**
	 * Gets the vms destroyed.
	 * 
	 * @return the vms destroyed
	 */
	protected int getVmsDestroyed() {
		return vmsDestroyed;
	}

	/**
	 * Sets the vms destroyed.
	 * 
	 * @param vmsDestroyed the new vms destroyed
	 */
	protected void setVmsDestroyed(int vmsDestroyed) {
		this.vmsDestroyed = vmsDestroyed;
	}

	/**
	 * Gets the datacenter ids list.
	 * 
	 * @return the datacenter ids list
	 */
	protected List<Integer> getDatacenterIdsList() {
		return datacenterIdsList;
	}

	/**
	 * Sets the datacenter ids list.
	 * 
	 * @param datacenterIdsList the new datacenter ids list
	 */
	protected void setDatacenterIdsList(List<Integer> datacenterIdsList) {
		this.datacenterIdsList = datacenterIdsList;
	}

	/**
	 * Gets the vms to datacenters map.
	 * 
	 * @return the vms to datacenters map
	 */
	protected Map<Integer, Integer> getVmsToDatacentersMap() {
		return vmsToDatacentersMap;
	}

	/**
	 * Sets the vms to datacenters map.
	 * 
	 * @param vmsToDatacentersMap the vms to datacenters map
	 */
	protected void setVmsToDatacentersMap(Map<Integer, Integer> vmsToDatacentersMap) {
		this.vmsToDatacentersMap = vmsToDatacentersMap;
	}

	/**
	 * Gets the datacenter characteristics list.
	 * 
	 * @return the datacenter characteristics list
	 */
	protected Map<Integer, DatacenterCharacteristics> getDatacenterCharacteristicsList() {
		return datacenterCharacteristicsList;
	}

	/**
	 * Sets the datacenter characteristics list.
	 * 
	 * @param datacenterCharacteristicsList the datacenter characteristics list
	 */
	protected void setDatacenterCharacteristicsList(
			Map<Integer, DatacenterCharacteristics> datacenterCharacteristicsList) {
		this.datacenterCharacteristicsList = datacenterCharacteristicsList;
	}

	/**
	 * Gets the datacenter requested ids list.
	 * 
	 * @return the datacenter requested ids list
	 */
	protected List<Integer> getDatacenterRequestedIdsList() {
		return datacenterRequestedIdsList;
	}

	/**
	 * Sets the datacenter requested ids list.
	 * 
	 * @param datacenterRequestedIdsList the new datacenter requested ids list
	 */
	protected void setDatacenterRequestedIdsList(List<Integer> datacenterRequestedIdsList) {
		this.datacenterRequestedIdsList = datacenterRequestedIdsList;
	}

}
