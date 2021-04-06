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
		//设定粗网格参数
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

		//生成3D网格
		GModel::current()->mesh(3);

		//区域分解
		GModel::current()->partitionMesh(
			numPart >= 0 ? numPart : CTX::instance()->mesh.numPartitions);

		//删除体网格
		std::for_each(GModel::current()->firstRegion(),
			GModel::current()->lastRegion(), deMeshGRegion());

		//加密
		for (int i = 0; i < refine; i++) {
			GModel::current()->refineMesh(CTX::instance()->mesh.secondOrderLinear,
				CTX::instance()->mesh.algoSubdivide == 2,
				CTX::instance()->mesh.algoSubdivide == 1,
				CTX::instance()->mesh.algoSubdivide == 3);
		}
		CTX::instance()->mesh.partitionSplitMeshFiles = 1;
		CTX::instance()->mesh.saveAll = 1;

		GmshWriteFile(partition_name + ".msh");//生成numpart份加密后的面网格
	}

	MPI_Barrier(MPI_COMM_WORLD);//同步语句，所有核等待上面主进程（0进程）生成numpart份.msh文件作为下面各个核的输入

	std::string str_id = std::to_string(id + 1);
	std::string str_temp = partition_name + "_" + str_id + ".msh";

	GmshOpenProject(str_temp);//每个核打开各自的.msh文件

	GModel::current()->mesh(3);//每个.msh文件都是加密后的面网格消息，根据面网格生成体网格

	CTX::instance()->mesh.partitionSplitMeshFiles = 0;

	GmshWriteFile("final_" + str_temp);

	std::cout << "------" + str_temp + "------" << std::endl;

	MPI_Barrier(MPI_COMM_WORLD);//同步，等待所有核生成好体网格

	//主进程合并各个体网格，生成最终文件final.msh
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