bool completed = false;
bool possible = true;
int numofvisited = 0;
int priorityQueue[M * N];

for (int z = 0; z < M * N; z++)
{
	priorityQueue[z] = -1;
}

//set all node distances to infinte and that they haven't been visited, also assigner there location with top 
//left being the origin

for (int i = 0; i < M; i++) {
	for (int j = 0; j < N; j++) {
		nodeArr[i * M + j].nodeEdit(i, j, INFINITE, false);
	}
}

//sets the start node within the array to have a distance of zero
nodeArr[start.getXCoord() * M + start.getYCoord()].setDis(0);

//This array will hold the order of nodes that need to be visited, zero means ignore that element
priorityQueue[0] = start.getXCoord() * M + start.getYCoord();

//check neighbouring elements to see if they're empty
int xVal;
int yVal;

xVal = nodeArr[priorityQueue[0]].getYCoord();
yVal = nodeArr[priorityQueue[0]].getXCoord();

while (completed == false || possible == true)
{
	//check if we are at the upper boarder
	if (yVal != 0)
	{
		printf("check above \n");
		if (maze[yVal - 1][xVal] == 'e')
		{
			nodeArr[(yVal - 1) * M + xVal].setDis(nodeArr[xVal * M + yVal].getDis() + 1);

			bool comp = false;
			int count = 0;
			int compNode = xVal * M + yVal;
			printf("above value node %d \n", compNode);
			int temp;

			//compare that the distance with the others and organise nodes accordingly
			while (comp == false)
			{
				if (nodeArr[priorityQueue[count]].getDis() > nodeArr[compNode].getDis())
				{
					temp = priorityQueue[count];
					priorityQueue[count] = compNode;
					compNode = temp;
				}
				else
				{
					comp = true;
				}
				count++;
			}
			if (maze[yVal - 1][xVal] == 'g') completed = true;
		}
	}

	if (yVal != M - 1)
	{
		printf("check below\n");
		if (maze[yVal + 1][xVal] == 'e')
		{
			nodeArr[(xVal + 1) * M + yVal].setDis(nodeArr[xVal * M + yVal].getDis() + 1);

			bool comp = false;
			int count = 0;
			int compNode = xVal * M + yVal;
			printf("below value node %d \n", compNode);
			int temp;

			//compare that the distance with the others
			while (comp == false)
			{
				if (nodeArr[priorityQueue[count]].getDis() > nodeArr[compNode].getDis())
				{
					temp = priorityQueue[count];
					priorityQueue[count] = compNode;
					compNode = temp;
				}
				else
				{
					comp = true;
				}
				count++;
			}
			if (maze[yVal + 1][xVal] == 'g') completed = true;
		}
	}

	if (xVal != 0)
	{
		printf("check right\n");
		if (maze[yVal][xVal - 1] == 'e')
		{
			nodeArr[(xVal - 1) * M + yVal].setDis(nodeArr[xVal * M + yVal].getDis() + 1);

			bool comp = false;
			int count = 0;
			int compNode = xVal * M + yVal;
			printf("right value node %d \n", compNode);
			int temp;

			//compare that the distance with the others
			while (comp == false)
			{
				if (nodeArr[priorityQueue[count]].getDis() > nodeArr[compNode].getDis())
				{
					temp = priorityQueue[count];
					priorityQueue[count] = compNode;
					compNode = temp;
				}
				else
				{
					comp = true;
				}
				count++;
			}
			if (maze[yVal][xVal - 1] == 'g') completed = true;
		}
	}

	if (xVal != N - 1)
	{
		printf("check left \n");
		if (maze[yVal][xVal + 1] == 'e')
		{
			nodeArr[(xVal + 1) * M + yVal].setDis(nodeArr[xVal * M + yVal].getDis() + 1);

			bool comp = false;
			int count = 0;
			int compNode = xVal * M + yVal;
			printf("left value node %d \n", compNode);
			int temp;

			//compare that the distance with the others
			while (comp == false)
			{
				if (nodeArr[priorityQueue[count]].getDis() > nodeArr[compNode].getDis())
				{
					temp = priorityQueue[count];
					priorityQueue[count] = compNode;
					compNode = temp;
				}
				else
				{
					comp = true;
				}
				count++;
			}
			if (maze[yVal][xVal + 1] == 'g') completed = true;
		}
	}
	xVal = nodeArr[priorityQueue[0]].getXCoord();
	yVal = nodeArr[priorityQueue[0]].getYCoord();

	for (int z = 0; z < M * N; z++)
	{
		printf("%d \n", priorityQueue[z]);
	}

}

/*class node {

	private:
		int x;
		int y;
		int distance;
		bool visited;
		//conNode connected;

	public:

		node() 
		{
			x = INFINITE;
			y = INFINITE;
			distance = INFINITE;
			visited = false;
		}

		node(int xVal, int yVal) 
		{
			x = xVal;
			y = yVal;
			distance = INFINITE;
			visited = false;
		}

		node(int xVal, int yVal, int dis, bool vis)
		{
			x = xVal;
			y = yVal;
			distance = dis;
			visited = vis;
		}
		
		int getXCoord()
		{
			return x;
		}

		int getYCoord()
		{
			return y;
		}

		int getDis() 
		{
			return distance;
		}

		bool getVis() 
		{
			return visited;
		}

		void setXVal(int xVal) 
		{
			x = xVal;
		}

		void setYVal(int yVal)
		{
			y = yVal;
		}

		void setDis(int dis)
		{
			distance = dis;
		}

		void setVis(bool vis)
		{
			visited = vis;
		}

		void nodeEdit(int xVal, int yVal)
		{
			x = xVal;
			y = yVal;
			distance = INFINITE;
			visited = false;
		}


		void nodeEdit(int xVal, int yVal, int dis, bool vis ) 
		{
			x = xVal;
			y = yVal;
			distance = dis;
			visited = vis;
		}
};*/