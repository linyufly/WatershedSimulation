/*
Program		:	Vincent-Soille Watershed Algorithm
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
/*
	printf("size of double *** = %lu\n", sizeof(double ***));

	printf("cube addr = %lu\n", cube);
	printf("&cube[0] = %lu\n", &cube[0]);
	printf("cube[0] addr = %lu\n", cube[0]);
	printf("cube[0][0] addr = %lu\n", cube[0][0]);
	printf("cube[0][0][0] val = %lf\n", cube[0][0][0]);
	printf("&cube[0][0][0] = %lu\n", &cube[0][0][0]);
*/
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

const int laplacianDelta = 16;

double ***ftleValues;
double spacing[3];
double origin[3];
int dimensions[3];
double hMax, hMin;

struct HeightComparator {
	bool operator () (const CellIndex &a, const CellIndex &b) const {
		return ftleValues[a.x][a.y][a.z] < ftleValues[b.x][b.y][b.z];
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

#define INIT -1
#define MASK -2
#define WSHED 0
#define FICTITIOUS CellIndex(-1, -1, -1)

void VincentSoille() {
	CellIndex *cellOrder = new CellIndex[dimensions[0] * dimensions[1] * dimensions[2]];
	int cnt = 0;
	for (int i = 0; i < dimensions[0]; i++)
		for (int j = 0; j < dimensions[1]; j++)
			for (int k = 0; k < dimensions[2]; k++)
				cellOrder[cnt++] = CellIndex(i, j, k);
	std::sort(cellOrder, cellOrder + cnt, HeightComparator());

	int ***lab = CreateCube<int>(dimensions[0], dimensions[1], dimensions[2]);
	int ***dist = CreateCube<int>(dimensions[0], dimensions[1], dimensions[2]);

	std::queue<CellIndex> queue;

	int curlab = 0;

	for (int i = 0; i < dimensions[0]; i++)
		for (int j = 0; j < dimensions[1]; j++)
			for (int k = 0; k < dimensions[2]; k++) {
				lab[i][j][k] = INIT;
				dist[i][j][k] = 0;
			}

	for (int d = 0; d < dimensions[0] * dimensions[1] * dimensions[2]; ) {
		//if (d >= 8000) break; // 1000, 4000, 6000

		CellIndex idx = cellOrder[d];
		double height = ftleValues[idx.x][idx.y][idx.z];

		int nextD;
		for (nextD = d; nextD < dimensions[0] * dimensions[1] * dimensions[2] &&
				ftleValues[cellOrder[nextD].x][cellOrder[nextD].y][cellOrder[nextD].z] == height;
				nextD++);
		
		for (int i = d; i < nextD; i++) {
			idx = cellOrder[i];
			lab[idx.x][idx.y][idx.z] = MASK;

			bool flag = false;
			for (int k = 0; k < 6; k++) {
				int _x = idx.x + dire[k][0];
				int _y = idx.y + dire[k][1];
				int _z = idx.z + dire[k][2];

				if (Outside(_x, _y, _z)) continue;

				if (lab[_x][_y][_z] > 0 || lab[_x][_y][_z] == WSHED) {
					flag = true;
					break;
				}
			}

			if (flag) {
				dist[idx.x][idx.y][idx.z] = 1;
				queue.push(idx);
			}
		}

		int curdist = 1;
		queue.push(FICTITIOUS);

		while (1) {
			idx = queue.front();
			queue.pop();

			if (idx.x == -1) {
				if (queue.empty()) break;
				else {
					queue.push(FICTITIOUS);
					curdist++;
				}
				idx = queue.front();
				queue.pop();
			}

			for (int k = 0; k < 6; k++) {
				int _x = idx.x + dire[k][0];
				int _y = idx.y + dire[k][1];
				int _z = idx.z + dire[k][2];

				if (Outside(_x, _y, _z)) continue;

				if (dist[_x][_y][_z] < curdist && (lab[_x][_y][_z] > 0 || lab[_x][_y][_z] == WSHED)) {
					/// DEBUG ///
					if (dist[_x][_y][_z] != curdist - 1) {
						printf("Found unexpected dist\n");
						exit(0);
					}

					if (lab[_x][_y][_z] > 0)
						if (lab[idx.x][idx.y][idx.z] == MASK || lab[idx.x][idx.y][idx.z] == WSHED)
							lab[idx.x][idx.y][idx.z] = lab[_x][_y][_z];
						else if (lab[idx.x][idx.y][idx.z] != lab[_x][_y][_z])
							lab[idx.x][idx.y][idx.z] = WSHED;
					else // lab[_x][_y][_z] == WSHED
						if (lab[idx.x][idx.y][idx.z] == MASK)
							lab[idx.x][idx.y][idx.z] = WSHED;
				} else if (lab[_x][_y][_z] == MASK && dist[_x][_y][_z] == 0) { // (_x, _y, _z) is plateau voxel
					dist[_x][_y][_z] = curdist + 1;
					queue.push(CellIndex(_x, _y, _z));
				}
			}
		}

		// New minima at level current height
		for (int i = d; i < nextD; i++) {
			idx = cellOrder[i];
			dist[idx.x][idx.y][idx.z] = 0;

			if (lab[idx.x][idx.y][idx.z] == MASK) {
				curlab++;
				queue.push(idx);
				lab[idx.x][idx.y][idx.z] = curlab;
				while (!queue.empty()) {
					idx = queue.front();
					queue.pop();
					for (int k = 0; k < 6; k++) {
						int _x = idx.x + dire[k][0];
						int _y = idx.y + dire[k][1];
						int _z = idx.z + dire[k][2];

						if (Outside(_x, _y, _z)) continue;

						if (lab[_x][_y][_z] == MASK) {
							queue.push(CellIndex(_x, _y, _z));
							lab[_x][_y][_z] = curlab;
						}
					}
				}
			}
		}

		d = nextD;
	}

	printf("curlab = %d\n", curlab);

	int numOfWatershedPixels = 0;

	for (int i = 0; i < dimensions[0]; i++)
		for (int j = 0; j < dimensions[1]; j++)
			for (int k = 0; k < dimensions[2]; k++)
				if (lab[i][j][k] < 0)
					printf("lab[%d][%d][%d] = %d\n", i, j, k, lab[i][j][k]);
				else if (lab[i][j][k] == 0) numOfWatershedPixels++;

	printf("numOfWatershedPixels = %d\n", numOfWatershedPixels);

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
				*pixel = lab[i][j][k];

				/// DEBUG ///
				if (lab[i][j][k] == WSHED || lab[i][j][k] == INIT) *pixel = curlab * 2;
			}

	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName("watershed_output.vtk");	

#if VTK_MAJOR_VERSION <= 5
	writer->SetInputConnection(image->GetProducerPort());
#else
	writer->SetInputData(image);
#endif

	writer->Write();

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

				/// DEBUG ///
				if (lab[i][j][k] == WSHED || lab[i][j][k] == INIT) *pixel = 0.0;
				//printf("%lf\n", *pixel);
			}

	writer->SetFileName("ftle.vtk");
	writer->Write();

	DeleteCube(lab);
	DeleteCube(dist);

	delete [] cellOrder;
}

int main() {
	/*
	double ***cube = CreateCube(100, 100, 100);
	DeleteCube(cube);
	*/
	LoadGrid();
	LaplacianSmoothing();
	VincentSoille();
	
	DeleteCube(ftleValues);

	return 0;
}
