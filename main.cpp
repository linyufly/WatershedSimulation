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

const int laplacianDelta = 0;

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
	reader->SetFileName("output_200.vtk");
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

void Simulation() {
	int ***labels = CreateCube<int>(dimensions[0], dimensions[1], dimensions[2]);
	for (int i = 0; i < dimensions[0]; i++)
		for (int j = 0; j < dimensions[1]; j++)
			for (int k = 0; k < dimensions[2]; k++)
				labels[i][j][k] = -1;

	priority_queue<CellIndex, std::vector<CellIndex>, HeightComparator> heap;

	for (int x = 0; x < dimensions[0]; x++)
		for (int y = 0; y < dimensions[1]; y++)
			for (int z = 0; z < dimensions[2]; z++) {
				if (labels[x][y][z] != -1) continue;
				bool flag = false;
				for (int k = 0; k < numOfNeighbors; k++) {
					int _x = x + dire[k][0];
					int _y = y + dire[k][1];
					int _z = z + dire[k][2];
					if (Outside(_x, _y, _z)) continue;
					if (ftleValues[_x][_y][_z] < ftleValues[x][y][z]) {
						flag = true;
						break;
					}
				}
				if (!flag) {
					FloodFill
				}
			}

	DeleteCube(labels);
}

int main() {
	LoadGrid();
	Simulation();
	
	DeleteCube(ftleValues);

	return 0;
}
