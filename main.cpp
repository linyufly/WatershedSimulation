/*
Program		:	Watershed Simulation
Author		:	Mingcheng Chen
*/

#include <vtkVersion.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>

#include <cstdio>

#include <queue>
#include <limits>
#include <algorithm>

template <class T>
T ***CreateCube(int X, int Y, int Z) {
	T *data = new T [X * Y * Z];
	T **zPointer = new T * [X * Y];
	for (int i = 0; i < X * Y; i++)
		zPointer[i] = data + Z * i;
	T ***yPointer = new T ** [X];
	for (int i = 0; i < X; i++)
		yPointer[i] = zPointer + Y * i;
	return yPointer;
}

template <class T>
void DeleteCube(T ***cube) {
	delete [] cube[0][0];
	delete [] cube[0];
	delete [] cube;
}

struct CellIndex {
	int x, y, z;

	CellIndex() {
	}

	CellIndex(int x, int y, int z) : x(x), y(y), z(z) {
	}
};

const int dire[6][3] = {{-1, 0, 0}, {1, 0, 0},
			{0, -1, 0}, {0, 1, 0},
			{0, 0, -1}, {0, 0, 1}};
const int numOfNeighbors = 6;

const int laplacianDelta = 4;

double ***ftleValues;
double spacing[3];
double origin[3];
int dimensions[3];
double hMax, hMin;

struct HeightComparator {
	bool operator () (const CellIndex &a, const CellIndex &b) const {
		return ftleValues[a.x][a.y][a.z] > ftleValues[b.x][b.y][b.z];
	}
};

bool Outside(int x, int y, int z) {
	return x < 0 || y < 0 || z < 0 || x >= dimensions[0] || y >= dimensions[1] || z >= dimensions[2];
}

void LoadGrid() {
	vtkSmartPointer<vtkStructuredPointsReader> reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
	reader->SetFileName("../VincentSoille/output.vtk");
	reader->Update();

	vtkStructuredPoints *structPoints = reader->GetOutput();
	structPoints->GetDimensions(dimensions);
	structPoints->GetOrigin(origin);
	structPoints->GetSpacing(spacing);

	printf("origin: %lf %lf %lf\n", origin[0], origin[1], origin[2]);
	printf("spacing: %lf %lf %lf\n", spacing[0], spacing[1], spacing[2]);

	ftleValues = CreateCube<double>(dimensions[0], dimensions[1], dimensions[2]);

	hMax = std::numeric_limits<double>::min();
	hMin = std::numeric_limits<double>::max();

	for (int i = 0; i < dimensions[0]; i++)
		for (int j = 0; j < dimensions[1]; j++)
			for (int k = 0; k < dimensions[2]; k++) {
				int pointId = i + j * dimensions[0] + k * dimensions[0] * dimensions[1];
				ftleValues[i][j][k] = structPoints->GetPointData()->GetScalars()->GetTuple1(pointId);
				if (ftleValues[i][j][k] > hMax) hMax = ftleValues[i][j][k];
				if (ftleValues[i][j][k] < hMin) hMin = ftleValues[i][j][k];
			}

	printf("hMax = %lf, hMin = %lf\n", hMax, hMin);

	/// DEBUG ///
	dimensions[2] = 1;
}

void LaplacianSmoothing() {
	printf("LaplacianSmoothing()\n");

	double ***tempValues = CreateCube<double>(dimensions[0], dimensions[1], dimensions[2]);

	for (int i = 0; i < dimensions[0]; i++)
		for (int j = 0; j < dimensions[1]; j++)
			for (int k = 0; k < dimensions[2]; k++) {
				tempValues[i][j][k] = 0;
				int cnt = 0;
				for (int dx = -laplacianDelta; dx <= laplacianDelta; dx++)
					for (int dy = -laplacianDelta; dy <= laplacianDelta; dy++)
						for (int dz = -laplacianDelta; dz <= laplacianDelta; dz++) {
							int x = i + dx;
							int y = j + dy;
							int z = k + dz;

							if (Outside(x, y, z)) continue;

							cnt++;
							tempValues[i][j][k] += ftleValues[x][y][z];
						}
				tempValues[i][j][k] /= cnt;
			}

	std::swap(ftleValues, tempValues);

	DeleteCube(tempValues);

	printf("Done.\n\n");
}

bool CheckMinima(int startX, int startY, int startZ, int nextLabel, int ***labels, bool ***visited, CellIndex *queue) {
	/// DEBUG ///
	//printf("(%d, %d, %d), nextLabel = %d\n", startX, startY, startZ, nextLabel);

	bool flag = false;
	int head, tail;
	queue[0] = CellIndex(startX, startY, startZ);
	visited[startX][startY][startZ] = true;
	for (head = tail = 0; head <= tail; head++) {
		/// DEBUG ///
		//printf("head = %d, tail = %d\n", head, tail);

		CellIndex currCell = queue[head];
		for (int k = 0; k < numOfNeighbors; k++) {
			int _x = currCell.x + dire[k][0];
			int _y = currCell.y + dire[k][1];
			int _z = currCell.z + dire[k][2];
			if (Outside(_x, _y, _z)) continue;
			if (ftleValues[_x][_y][_z] < ftleValues[startX][startY][startZ]) flag = true;

			if (visited[_x][_y][_z]) continue;
			if (ftleValues[_x][_y][_z] == ftleValues[startX][startY][startZ]) {
				visited[_x][_y][_z] = true;
				queue[++tail] = CellIndex(_x, _y, _z);
			}
		}
	}
	if (!flag) {
		for (int i = 0; i <= tail; i++)
			labels[queue[i].x][queue[i].y][queue[i].z] = nextLabel;
		return true;
	}
	return false;
}

void Simulation() {
	printf("Simulation()\n");

	int ***labels = CreateCube<int>(dimensions[0], dimensions[1], dimensions[2]);
	for (int i = 0; i < dimensions[0]; i++)
		for (int j = 0; j < dimensions[1]; j++)
			for (int k = 0; k < dimensions[2]; k++)
				labels[i][j][k] = -1;

	bool ***visited = CreateCube<bool>(dimensions[0], dimensions[1], dimensions[2]);
	for (int i = 0; i < dimensions[0]; i++)
		for (int j = 0; j < dimensions[1]; j++)
			for (int k = 0; k < dimensions[2]; k++)
				visited[i][j][k] = false;

	/// DEBUG ///
	printf("before finding minima\n");

	int nextLabel = 0;

	CellIndex *queueArray = new CellIndex [dimensions[0] * dimensions[1] * dimensions[2]];

	for (int x = 0; x < dimensions[0]; x++)
		for (int y = 0; y < dimensions[1]; y++)
			for (int z = 0; z < dimensions[2]; z++)
				if (!visited[x][y][z])
					if (CheckMinima(x, y, z, nextLabel, labels, visited, queueArray)) nextLabel++;

	delete [] queueArray;

	/// DEBUG ///
	printf("nextLabel = %d\n", nextLabel);
	printf("Finish finding minima\n");

	std::priority_queue<CellIndex, std::vector<CellIndex>, HeightComparator> heap;

	for (int i = 0; i < dimensions[0]; i++)
		for (int j = 0; j < dimensions[1]; j++)
			for (int k = 0; k < dimensions[2]; k++)
				if (labels[i][j][k] != -1) {
					bool flag = false;
					for (int d = 0; d < numOfNeighbors; d++) {
						int _x = i + dire[d][0];
						int _y = j + dire[d][1];
						int _z = k + dire[d][2];
						if (Outside(_x, _y, _z)) continue;
						if (labels[_x][_y][_z] == -1) {
							flag = true;
							break;
						}
					}
					if (flag) heap.push(CellIndex(i, j, k));
				}

	while (!heap.empty()) {
		CellIndex curr = heap.top();
		heap.pop();
		for (int k = 0; k < numOfNeighbors; k++) {
			int _x = curr.x + dire[k][0];
			int _y = curr.y + dire[k][1];
			int _z = curr.z + dire[k][2];
			if (Outside(_x, _y, _z)) continue;
			if (labels[_x][_y][_z] != -1) continue;
			if (ftleValues[_x][_y][_z] >= ftleValues[curr.x][curr.y][curr.z]) {
				labels[_x][_y][_z] = labels[curr.x][curr.y][curr.z];
				heap.push(CellIndex(_x, _y, _z));
			}
		}
	}

	int *labelColors = new int [nextLabel];
	for (int i = 0; i < nextLabel; i++)
		labelColors[i] = i;
	std::random_shuffle(labelColors, labelColors + nextLabel);

	vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
	image->SetDimensions(dimensions);
	image->SetOrigin(origin);
	image->SetSpacing(spacing);

#if VTK_MAJOR_VERSION <= 5
	image->SetNumberOfScalarComponents(1);
	image->SetScalarTypeToInt();
#else
	image->AllocateScalars(VTK_INT, 1);
#endif

	for (int i = 0; i < dimensions[0]; i++)
		for (int j = 0; j < dimensions[1]; j++)
			for (int k = 0; k < dimensions[2]; k++) {
				int *pixel = static_cast<int *>(image->GetScalarPointer(i, j, k));
				*pixel = labels[i][j][k];
				if (labels[i][j][k] >= 0) *pixel = labelColors[*pixel];
			}

	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName("watershed_output.vtk");	

#if VTK_MAJOR_VERSION <= 5
	writer->SetInputConnection(image->GetProducerPort());
#else
	writer->SetInputData(image);
#endif

	writer->Write();

	delete [] labelColors;

#if VTK_MAJOR_VERSION <= 5
	image->SetNumberOfScalarComponents(1);
	image->SetScalarTypeToDouble();
#else
	image->AllocateScalars(VTK_DOUBLE, 1);
#endif

	for (int i = 0; i < dimensions[0]; i++)
		for (int j = 0; j < dimensions[1]; j++)
			for (int k = 0; k < dimensions[2]; k++) {
				double *pixel = static_cast<double *>(image->GetScalarPointer(i, j, k));
				*pixel = ftleValues[i][j][k];
			}

	writer->SetFileName("ftle.vtk");
	writer->Write();

	DeleteCube(labels);
	DeleteCube(visited);

	printf("Done\n");
}

int main() {
	LoadGrid();
	LaplacianSmoothing();
	Simulation();
	
	DeleteCube(ftleValues);

	return 0;
}
