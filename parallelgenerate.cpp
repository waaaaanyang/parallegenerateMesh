//#if defined(HAVE_MPI)
GMSH_API void gmsh::model::mesh::parallelgenerate(
	const double size, const int numPart, const int refine,
	const std::string &fileName,
	const std::string &partition_name) 
{
	int id;
	int p;

	MPI_Init(&_argc, &_argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	if (id == 0) {
		GmshOpenProject(fileName + ".geo");
		//�趨���������
		std::vector<std::pair<int, int> > tmp;
		tmp.clear();
		std::vector<GEntity *> entities;
		GModel::current()->getEntities(entities, 0);
		for (std::size_t i = 0; i < entities.size(); i++)
			tmp.push_back(
				std::pair<int, int>(entities[i]->dim(), entities[i]->tag()));

		for (std::size_t i = 0; i < tmp.size(); i++) {
			int dim = tmp[i].first, tag = tmp[i].second;
			if (dim == 0) {
				GVertex *gv = GModel::current()->getVertexByTag(tag);
				if (gv) gv->setPrescribedMeshSizeAtVertex(size);
			}
		}

		//����3D����
		GModel::current()->mesh(3);

		//����ֽ�
		GModel::current()->partitionMesh(
			numPart >= 0 ? numPart : CTX::instance()->mesh.numPartitions);

		//ɾ��������
		std::for_each(GModel::current()->firstRegion(),
			GModel::current()->lastRegion(), deMeshGRegion());

		//����
		for (int i = 0; i < refine; i++) {
			GModel::current()->refineMesh(CTX::instance()->mesh.secondOrderLinear,
				CTX::instance()->mesh.algoSubdivide == 2,
				CTX::instance()->mesh.algoSubdivide == 1,
				CTX::instance()->mesh.algoSubdivide == 3);
		}
		CTX::instance()->mesh.partitionSplitMeshFiles = 1;
		CTX::instance()->mesh.saveAll = 1;

		GmshWriteFile(partition_name + ".msh");//����numpart�ݼ��ܺ��������
	}

	MPI_Barrier(MPI_COMM_WORLD);//ͬ����䣬���к˵ȴ����������̣�0���̣�����numpart��.msh�ļ���Ϊ��������˵�����

	std::string str_id = std::to_string(id + 1);
	std::string str_temp = partition_name + "_" + str_id + ".msh";

	GmshOpenProject(str_temp);//ÿ���˴򿪸��Ե�.msh�ļ�

	GModel::current()->mesh(3);//ÿ��.msh�ļ����Ǽ��ܺ����������Ϣ����������������������

	CTX::instance()->mesh.partitionSplitMeshFiles = 0;

	GmshWriteFile("final_" + str_temp);

	std::cout << "------" + str_temp + "------" << std::endl;

	MPI_Barrier(MPI_COMM_WORLD);//ͬ�����ȴ����к����ɺ�������

	//�����̺ϲ��������������������ļ�final.msh
	if (id == 0) {

		GmshOpenProject("final_" + partition_name + "_" + "1.msh");

		for (int i = 2; i <= numPart; i++) {
			std::string str_i = std::to_string(i);
			GmshMergeFile("final_" + partition_name + "_" + str_i + ".msh");
		}

		GModel::current()->unpartitionMesh();

		CTX::instance()->mesh.mshFileVersion = 2.0;

		GmshWriteFile("final.msh");

	}

	MPI_Finalize();
}
//#endif