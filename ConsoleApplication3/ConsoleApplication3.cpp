
#define OLC_PGE_APPLICATION

#include "olcPixelGameEngine.h"

#include <iostream>
#include <stack>
#include <string>
#include <algorithm>
#include <chrono>
#include <thread>

#define SCREENHEIGHT 720
#define SCREENWIDTH 720

#define MAZEWIDTH 100
#define MAZEHEIGHT 100

//primnodes mazeNodes;
olc::Pixel backColour = olc::BLACK;

// Override base class with your custom functionality
class pathFindingAlgorithm : public olc::PixelGameEngine
{

public:

	pathFindingAlgorithm()
	{
		// Name you application
		sAppName = "Dijkstras Algorithm";
	}

private:

	struct sCell					//used to generate maze
	{
		int x = 0;															//x position of cell or collumn value in 2D matrix
		int y = 0;															//y position of cell or row value in 2D matrix
		char cCellRep = 'w';												//the letter which tell us what it is
		bool bVisitedOne = false;											//This will let us know the first time the cell is visited
		bool bVisitedTwo = false;											//This will be used as an error checker as if a cell is visited more than twice an error has occured
		bool NORTH = false;
		bool SOUTH = false;
		bool EAST = false;
		bool WEST = false;
	};

	struct sNode					//used to solve maze
	{
		bool bObstacle = false;												//checks if this node is a wall
		bool bVisited = false;												//checks if node has been visited
		float fStartToGoal = 0.0f;											//summation of distance from start to node and heuristic (Added for A*)
		float fStartToNode = 0.0f;											//distance from start node to current node (Djikstra)
		int x = 0;															//x position
		int y = 0;															//y position
		std::vector<sNode*> vecNeighbours;									//vector containing pointer to the neighbouring nodes
		sNode* parent = nullptr;											//pointer to parent node(node which was last to update it)
	};

	sCell cellStart;														//starting cell
	sCell cellEnd;															//ending cell
	
	sCell* cells = nullptr;													//an array of cells that will be used to create our 2D array

	std::stack<sCell> recursiveStack;										//Create the stack for a recursive effect

	char maze[MAZEHEIGHT][MAZEWIDTH];										//This will be a 2D array that holds the generated map
	
	int numOfVisitedCells = 0;												//This will be used to cancel the loop if all cells within the array have been visited
	
	sNode* nodes = nullptr;													//initialise a pointer to an array of nodes	

	bool bRefreshPath = false;

	sNode* nodeStart = nullptr;												//starting node
	sNode* nodeEnd = nullptr;												//ending node

	int boxHeight = SCREENHEIGHT / MAZEHEIGHT;								//ratio which converts virtual to real screen
	int boxWidth = SCREENWIDTH / MAZEWIDTH;									//ratio which converts virtual to real screen

public:

	void createMaze()
	{
		//printf("\n inside the create maze function");

		//for loop to initialise 
		for (int i = 0; i < MAZEWIDTH; i++) {
			for (int j = 0; j < MAZEHEIGHT; j++)
			{
				//set cell coordinates
				cells[j * MAZEWIDTH + i].NORTH = false;
				cells[j * MAZEWIDTH + i].SOUTH = false;
				cells[j * MAZEWIDTH + i].EAST = false;
				cells[j * MAZEWIDTH + i].WEST = false;
				cells[j * MAZEWIDTH + i].bVisitedOne = false;
				cells[j * MAZEWIDTH + i].bVisitedTwo = false;
				cells[j * MAZEWIDTH + i].cCellRep = 'w';
				cells[j * MAZEWIDTH + i].x = i;
				cells[j * MAZEWIDTH + i].y = j;
			}
		}

		// randomly select starting node
		cellStart.x = rand() % MAZEWIDTH;
		cellStart.y = rand() % MAZEHEIGHT;
		cells[cellStart.y * MAZEWIDTH + cellStart.x].cCellRep = 's';

		//randomly select ending node
		cellEnd.x = rand() % MAZEWIDTH;
		cellEnd.y = rand() % MAZEHEIGHT;
		cells[cellEnd.y * MAZEWIDTH + cellEnd.x].cCellRep = 'g';

		//push start cell onto it
		recursiveStack.push(cells[cellStart.y * MAZEWIDTH + cellStart.x]);

		numOfVisitedCells = 0;

		int numOfNextCell = 0;													//This will hold a value 0-3 which corresponds to NORTH,EAST,SOUTH and WEST respectively
		int numOfCells = MAZEHEIGHT*MAZEWIDTH;

		while ((numOfVisitedCells < numOfCells) && (!recursiveStack.empty()))
		{

			//NORTH:0 EAST:1 SOUTH:2 WEST:3 random decision
			numOfNextCell = rand() % 4;

			switch (numOfNextCell) 
			{
			case 0: //NORTH
				cells[(recursiveStack.top().y) * MAZEWIDTH + recursiveStack.top().x].NORTH = true;
				break;
			case 1://EAST
				cells[(recursiveStack.top().y) * MAZEWIDTH + recursiveStack.top().x].EAST = true;
				break;
			case 2://SOUTH
				cells[(recursiveStack.top().y) * MAZEWIDTH + recursiveStack.top().x].SOUTH = true;
				break;
			case 3://WEST
				cells[(recursiveStack.top().y) * MAZEWIDTH + recursiveStack.top().x].WEST = true;
			}

			//printf("\n number of cells visited = %d", numOfVisitedCells);
			//printf("\n cell %d visited", recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x);
			//printf("\n random value = %d", numOfNextCell);
			//printf("\n number of elements in stack = %d", (int)recursiveStack.size());

			if (cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x].bVisitedOne == false)
			{
				cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x].bVisitedOne = true;
			}
			else
			{
				cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x].bVisitedTwo = true;
			}

			//select top cell
			if (recursiveStack.top().y > 0 &&
				numOfNextCell == 0 &&
				cells[(recursiveStack.top().y - 1) * MAZEWIDTH + recursiveStack.top().x].cCellRep == 'w' && cells[(recursiveStack.top().y - 1) * MAZEWIDTH + recursiveStack.top().x].cCellRep != 's' && //this checks the NORTH cell
				cells[(recursiveStack.top().y - 2) * MAZEWIDTH + recursiveStack.top().x].cCellRep == 'w' && cells[(recursiveStack.top().y - 2) * MAZEWIDTH + recursiveStack.top().x].cCellRep != 's' && //this checks the cell above the NORTH cell
				cells[(recursiveStack.top().y - 1) * MAZEWIDTH + recursiveStack.top().x - 1].cCellRep == 'w' && cells[(recursiveStack.top().y - 1) * MAZEWIDTH + recursiveStack.top().x - 1].cCellRep != 's' && //this checks the NW cell
				cells[(recursiveStack.top().y - 1) * MAZEWIDTH + recursiveStack.top().x + 1].cCellRep == 'w' && cells[(recursiveStack.top().y - 1) * MAZEWIDTH + recursiveStack.top().x + 1].cCellRep != 's' //this checks the NE cell
				//&& cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x].bVisitedTwo == false
				) {

				cells[(recursiveStack.top().y - 1) * MAZEWIDTH + recursiveStack.top().x].cCellRep = 'e';
				recursiveStack.push(cells[(recursiveStack.top().y - 1) * MAZEWIDTH + recursiveStack.top().x]);
				numOfVisitedCells++;
			}
			//select south cell
			else if (recursiveStack.top().y < MAZEHEIGHT - 1 &&
				numOfNextCell == 2 &&
				cells[(recursiveStack.top().y + 1) * MAZEWIDTH + recursiveStack.top().x].cCellRep == 'w' && cells[(recursiveStack.top().y + 1) * MAZEWIDTH + recursiveStack.top().x].cCellRep != 's' &&				//this checks the SOUTH cell
				cells[(recursiveStack.top().y + 2) * MAZEWIDTH + recursiveStack.top().x].cCellRep == 'w' && cells[(recursiveStack.top().y + 2) * MAZEWIDTH + recursiveStack.top().x].cCellRep != 's' &&				//this checks the cell below the SOUTH cell
				cells[(recursiveStack.top().y + 1) * MAZEWIDTH + recursiveStack.top().x - 1].cCellRep == 'w' && cells[(recursiveStack.top().y + 1) * MAZEWIDTH + recursiveStack.top().x - 1].cCellRep != 's' &&		//this checks the NW cell
				cells[(recursiveStack.top().y + 1) * MAZEWIDTH + recursiveStack.top().x + 1].cCellRep == 'w' && cells[(recursiveStack.top().y + 1) * MAZEWIDTH + recursiveStack.top().x + 1].cCellRep != 's' 		//this checks the NE cell
				//&& cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x].bVisitedTwo == false
				) {

				cells[(recursiveStack.top().y + 1) * MAZEWIDTH + recursiveStack.top().x].cCellRep = 'e';
				recursiveStack.push(cells[(recursiveStack.top().y + 1) * MAZEWIDTH + recursiveStack.top().x]);
				numOfVisitedCells++;
			}
			//select west cell
			else if (recursiveStack.top().x > 0 && numOfNextCell == 3 &&
				cells[(recursiveStack.top().y) * MAZEWIDTH + recursiveStack.top().x - 1].cCellRep == 'w' && cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x - 1].cCellRep != 's' &&				//This checks the cell to the EAST
				cells[(recursiveStack.top().y) * MAZEWIDTH + recursiveStack.top().x - 2].cCellRep == 'w' && cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x - 2].cCellRep != 's' &&				//This checks the cell to the EAST of the EAST cell
				cells[(recursiveStack.top().y - 1) * MAZEWIDTH + recursiveStack.top().x - 1].cCellRep == 'w' && cells[(recursiveStack.top().y - 1) * MAZEWIDTH + recursiveStack.top().x - 1].cCellRep != 's' &&		//this checks the NW cell
				cells[(recursiveStack.top().y + 1) * MAZEWIDTH + recursiveStack.top().x - 1].cCellRep == 'w' && cells[(recursiveStack.top().y + 1) * MAZEWIDTH + recursiveStack.top().x - 1].cCellRep != 's' 		//this checks the SW cell
				//&& cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x].bVisitedTwo == false
				) {

				cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x - 1].cCellRep = 'e';
				recursiveStack.push(cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x - 1]);
				numOfVisitedCells++;

			}
			//select east cell
			else if (recursiveStack.top().x < MAZEWIDTH - 1 && numOfNextCell == 1 &&
				cells[(recursiveStack.top().y) * MAZEWIDTH + recursiveStack.top().x + 1].cCellRep == 'w' && cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x + 1].cCellRep != 's' &&				//This checks the cell to the EAST
				cells[(recursiveStack.top().y) * MAZEWIDTH + recursiveStack.top().x + 2].cCellRep == 'w' && cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x + 2].cCellRep != 's' &&				//This checks the cell to the EAST of the EAST cell
				cells[(recursiveStack.top().y - 1) * MAZEWIDTH + recursiveStack.top().x + 1].cCellRep == 'w' && cells[(recursiveStack.top().y - 1) * MAZEWIDTH + recursiveStack.top().x + 1].cCellRep != 's' &&		//this checks the NW cell
				cells[(recursiveStack.top().y + 1) * MAZEWIDTH + recursiveStack.top().x + 1].cCellRep == 'w' && cells[(recursiveStack.top().y + 1) * MAZEWIDTH + recursiveStack.top().x + 1].cCellRep != 's' 		//this checks the SW cell
				//&& cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x].bVisitedTwo == false
				) {

				cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x + 1].cCellRep = 'e';
				recursiveStack.push(cells[recursiveStack.top().y * MAZEWIDTH + recursiveStack.top().x + 1]);
				numOfVisitedCells++;
			}
			else if (cells[(recursiveStack.top().y) * MAZEWIDTH + recursiveStack.top().x].NORTH == true &&
				cells[(recursiveStack.top().y) * MAZEWIDTH + recursiveStack.top().x].SOUTH == true &&
				cells[(recursiveStack.top().y) * MAZEWIDTH + recursiveStack.top().x].WEST == true &&
				cells[(recursiveStack.top().y) * MAZEWIDTH + recursiveStack.top().x].EAST == true) recursiveStack.pop();
			//if nothing is selected we have reached a dead end so we need to pop that cell and back track	
		}

		//printf("\n The stack is done");

		for (int i = 0; i < MAZEWIDTH; i++) {
			for (int j = 0; j < MAZEHEIGHT; j++) 
			{
				maze[j][i] = cells[j * MAZEWIDTH + i].cCellRep;
			}
		}

		bool bBlocked = true;

		while (bBlocked) {
			numOfNextCell = rand() % 4;

			switch (numOfNextCell)
			{
			case 0://NORTH
				if (cellEnd.y > 0) {
					if (maze[cellEnd.y-1][cellEnd.x] == 'w') {
						maze[cellEnd.y-1][cellEnd.x] = 'e';
						bBlocked = false;
					}
				}
				break;
			case 1://EAST
				if (cellEnd.x < MAZEWIDTH - 1) {
					if (maze[cellEnd.y][cellEnd.x + 1] == 'w') {
						maze[cellEnd.y][cellEnd.x + 1] = 'e';
						bBlocked = false;
					}
				}
				break;
			case 2://SOUTH
				if (cellEnd.y < MAZEHEIGHT -1) {
					if (maze[cellEnd.y+1][cellEnd.x] == 'w') {
						maze[cellEnd.y+1][cellEnd.x] = 'e';
						bBlocked = false;
					}
				}
				break;
			case 3://WEST
				if (cellEnd.x > 0) {
					if (maze[cellEnd.y][cellEnd.x-1] == 'w') {
						maze[cellEnd.y][cellEnd.x-1] = 'e';
						bBlocked = false;
					}
				}
				break;
			}
		}

		bBlocked = true;

		while (bBlocked) {
			numOfNextCell = rand() % 4;

			switch (numOfNextCell)
			{
			case 0://NORTH
				if (cellStart.y > 0) {
					if (maze[cellStart.y - 1][cellStart.x] == 'w') {
						maze[cellStart.y - 1][cellStart.x] = 'e';
						bBlocked = false;
					}
				}
				break;
			case 1://EAST
				if (cellStart.x < MAZEWIDTH - 1) {
					if (maze[cellStart.y][cellStart.x + 1] == 'w') {
						maze[cellStart.y][cellStart.x + 1] = 'e';
						bBlocked = false;
					}
				}
				break;
			case 2://SOUTH
				if (cellStart.y < MAZEHEIGHT - 1) {
					if (maze[cellStart.y + 1][cellStart.x] == 'w') {
						maze[cellStart.y + 1][cellStart.x] = 'e';
						bBlocked = false;
					}
				}
				break;
			case 3://WEST
				if (cellStart.x > 0) {
					if (maze[cellStart.y][cellStart.x - 1] == 'w') {
						maze[cellStart.y][cellStart.x - 1] = 'e';
						bBlocked = false;
					}
				}
				break;
			}
		}
	}

	bool solveDijkstra()
	{
		printf("\n solve via Djikstra");

		for (int i = 0; i < MAZEWIDTH; i++) {
			for (int j = 0; j < MAZEHEIGHT; j++)
			{
				nodes[j * MAZEWIDTH + i].bObstacle = false;
				nodes[j * MAZEWIDTH + i].bVisited = false;
				nodes[j * MAZEWIDTH + i].x = i;
				nodes[j * MAZEWIDTH + i].y = j;
				nodes[j * MAZEWIDTH + i].parent = nullptr;
			}
		}

		//establish connections between adjecents nodes
		for (int i = 0; i < MAZEWIDTH; i++) {
			for (int j = 0; j < MAZEHEIGHT; j++)
			{
				//if not along the outer border we will add the memory address of the cell that was just checked to the vector 
				if (j > 0) nodes[j * MAZEWIDTH + i].vecNeighbours.push_back(&nodes[(j - 1) * MAZEWIDTH + i + 0]);
				if (j < MAZEHEIGHT - 1) nodes[j * MAZEWIDTH + i].vecNeighbours.push_back(&nodes[(j + 1) * MAZEWIDTH + i + 0]);
				if (i > 0) nodes[j * MAZEWIDTH + i].vecNeighbours.push_back(&nodes[j * MAZEWIDTH + i - 1]);
				if (i < MAZEWIDTH - 1) nodes[j * MAZEWIDTH + i].vecNeighbours.push_back(&nodes[j * MAZEWIDTH + i + 1]);
			}
		}

		for (int i = 0; i < MAZEWIDTH; i++) {
			for (int j = 0; j < MAZEHEIGHT; j++)
			{
				if (maze[j][i] == 's') nodeStart = &nodes[j * MAZEWIDTH + i];
				if (maze[j][i] == 'g') nodeEnd = &nodes[j * MAZEWIDTH + i];
				if (maze[j][i] == 'w') nodes[j * MAZEWIDTH + i].bObstacle = true;
			}
		}

		//reset the nodes
		for (int i = 0; i < MAZEWIDTH; i++) {
			for (int j = 0; j < MAZEHEIGHT; j++)
			{
				nodes[j * MAZEWIDTH + i].bVisited = false;
				nodes[j * MAZEWIDTH + i].fStartToGoal = INFINITE;
				nodes[j * MAZEWIDTH + i].fStartToNode = INFINITE;
				nodes[j * MAZEWIDTH + i].parent = nullptr;
			}
		}

		sNode* nodeCurrent = nodeStart;
		nodeStart->fStartToNode = 0.0f;

		//this will use a heuristic for A*
		nodeEnd->fStartToGoal = 0.0f;

		std::list<sNode*> listUnvisitedNodes;
		listUnvisitedNodes.push_back(nodeStart);

		while (!listUnvisitedNodes.empty()) //&& nodeCurrent != nodeEnd
		{
			//sort the list by value
			listUnvisitedNodes.sort([](const sNode* a, const sNode* b) {return a->fStartToNode < b->fStartToNode; });

			//This will check if the any nodes have been visited and if so to remove them
			while (!listUnvisitedNodes.empty() && listUnvisitedNodes.front()->bVisited)
				listUnvisitedNodes.pop_front();

			//if the list is empty we want to stop the operation
			if (listUnvisitedNodes.empty())
				break;

			//sets node at the front of the list to be the next explored
			nodeCurrent = listUnvisitedNodes.front();
			//as it has now been visited we set it to true
			nodeCurrent->bVisited = true;

			//check
			for (auto nodeNeighbour : nodeCurrent->vecNeighbours)
			{
				if (!nodeNeighbour->bVisited && nodeNeighbour->bObstacle == false)
					listUnvisitedNodes.push_back(nodeNeighbour);

				float fDistanceCheck = nodeCurrent->fStartToNode;

				if (fDistanceCheck < nodeNeighbour->fStartToNode)
				{
					nodeNeighbour->parent = nodeCurrent;
					nodeNeighbour->fStartToNode = fDistanceCheck;
				}
			}

			//std::this_thread::sleep_for(std::chrono::milliseconds(10));

		}

		printf("\n solved");

		return true;
	}

	bool solveAStar()
	{
		printf("\n Start A*");

		for (int i = 0; i < MAZEWIDTH; i++) {
			for (int j = 0; j < MAZEHEIGHT; j++)
			{
				nodes[j * MAZEWIDTH + i].bObstacle = false;
				nodes[j * MAZEWIDTH + i].bVisited = false;
				nodes[j * MAZEWIDTH + i].x = i;
				nodes[j * MAZEWIDTH + i].y = j;
				nodes[j * MAZEWIDTH + i].parent = nullptr;
			}
		}

		//establish connections between adjecents nodes
		for (int i = 0; i < MAZEWIDTH; i++) {
			for (int j = 0; j < MAZEHEIGHT; j++)
			{
				//if not along the outer border we will add the memory address of the cell that was just checked to the vector 
				if (j > 0) nodes[j * MAZEWIDTH + i].vecNeighbours.push_back(&nodes[(j - 1) * MAZEWIDTH + i + 0]);
				if (j < MAZEHEIGHT - 1) nodes[j * MAZEWIDTH + i].vecNeighbours.push_back(&nodes[(j + 1) * MAZEWIDTH + i + 0]);
				if (i > 0) nodes[j * MAZEWIDTH + i].vecNeighbours.push_back(&nodes[j * MAZEWIDTH + i - 1]);
				if (i < MAZEWIDTH - 1) nodes[j * MAZEWIDTH + i].vecNeighbours.push_back(&nodes[j * MAZEWIDTH + i + 1]);
			}
		}

		for (int i = 0; i < MAZEWIDTH; i++) {
			for (int j = 0; j < MAZEHEIGHT; j++)
			{
				if (maze[j][i] == 's') nodeStart = &nodes[j * MAZEWIDTH + i];
				if (maze[j][i] == 'g') nodeEnd = &nodes[j * MAZEWIDTH + i];
				if (maze[j][i] == 'w') nodes[j * MAZEWIDTH + i].bObstacle = true;
			}
		}

		//reset the nodes
		for (int i = 0; i < MAZEWIDTH; i++) {
			for (int j = 0; j < MAZEHEIGHT; j++)
			{
				nodes[j * MAZEWIDTH + i].bVisited = false;
				nodes[j * MAZEWIDTH + i].fStartToGoal = INFINITE;
				nodes[j * MAZEWIDTH + i].fStartToNode = INFINITE;
				nodes[j * MAZEWIDTH + i].parent = nullptr;
			}
		}

		//Pythagoras to solve for euclidean displacement
		auto distance = [](sNode* a, sNode* b)
		{
			return sqrtf((a->x - b->x) * (a->x - b->x) + (a->y - b->y) * (a->y - b->y));
		};

		// for simplicity the heuristic will just use euclidean displacement
		auto heuristic = [distance](sNode* a, sNode* b)
		{
			return distance(a, b);
		};

		sNode* nodeCurrent = nodeStart;
		nodeStart->fStartToNode = 0.0f;
		nodeEnd->fStartToGoal = heuristic(nodeStart, nodeEnd);

		//creates a list of nodes and adds the start node to the top
		std::list<sNode*> listUnvisitedNodes;
		listUnvisitedNodes.push_back(nodeStart);

		while (!listUnvisitedNodes.empty() && nodeCurrent != nodeEnd)
		{
			//sort the list by value
			listUnvisitedNodes.sort([](const sNode* lhs, const sNode* rhs) {return lhs->fStartToGoal < rhs->fStartToGoal; });

			//This will check if the any nodes have been visited and if so to remove them
			while (!listUnvisitedNodes.empty() && listUnvisitedNodes.front()->bVisited)
				listUnvisitedNodes.pop_front();

			//if the list is empty we want to stop the operation
			if (listUnvisitedNodes.empty())
				break;

			//sets node at the front of the list to be the next explored
			nodeCurrent = listUnvisitedNodes.front();
			//as it has now been visited we set it to true
			nodeCurrent->bVisited = true;

			//for loop which checks each neighbour of the current node
			for (auto nodeNeighbour : nodeCurrent->vecNeighbours)
			{
				//if the neighbour is unvisited and not an obstacle will it be added to the list
				if (!nodeNeighbour->bVisited && nodeNeighbour->bObstacle == false)
					listUnvisitedNodes.push_back(nodeNeighbour);

				//this variable is used to compare the new distance from start to neighbour node with the neighbours previous distance
				float fDistanceCheck = nodeCurrent->fStartToNode + distance(nodeCurrent, nodeNeighbour);

				//if the distance is less we will enter this loop
				if (fDistanceCheck < nodeNeighbour->fStartToNode)
				{
					//this will set the parent of the neighbour node to the current as it provides a shorter distance
					nodeNeighbour->parent = nodeCurrent;
					//this will then set the distance from start to the neighbour node
					nodeNeighbour->fStartToNode = fDistanceCheck;
					//this will then add the heuristic for the total estimated distance via the neighbour route
					nodeNeighbour->fStartToGoal = nodeNeighbour->fStartToNode + heuristic(nodeNeighbour, nodeEnd);
				}
			}

			//std::this_thread::sleep_for(std::chrono::milliseconds(100));

		}

		printf("\n solved");

		return true;
	}

	bool OnUserCreate() override	// Called once at the start
	{
		cells = new sCell[MAZEHEIGHT * MAZEWIDTH];

		printf("\n Maze creation start");

		createMaze();

		printf("\n Maze created");
		//This initialises every node with their respective coordinates and sets their boolean variables
		nodes = new sNode[MAZEHEIGHT * MAZEWIDTH];

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{

		for (int i = 0; i < MAZEWIDTH; i++) {
			for (int j = 0; j < MAZEHEIGHT; j++)
			{
				int topLeft_Y = j * boxHeight;
				int topLeft_X = i * boxWidth;

				if (maze[j][i] == 'w')
				{
					FillRect(topLeft_X, topLeft_Y, boxWidth, boxHeight, backColour);
					FillRect(topLeft_X + 1, topLeft_Y + 1, boxWidth - 2, boxHeight - 2, olc::BLUE);
				}
				else if (maze[j][i] == 'e')
				{
					FillRect(topLeft_X, topLeft_Y, boxWidth, boxHeight, backColour);
					FillRect(topLeft_X + 1, topLeft_Y + 1, boxWidth - 2, boxHeight - 2, nodes[j * MAZEWIDTH + i].bVisited && bRefreshPath == true ? olc::GREY : olc::WHITE);
				}
				else if (maze[j][i] == 's')
				{
					FillRect(topLeft_X, topLeft_Y, boxWidth, boxHeight, backColour);
					FillRect(topLeft_X + 1, topLeft_Y + 1, boxWidth - 2, boxHeight - 2, olc::DARK_MAGENTA);
				}
				else if (maze[j][i] == 'g')
				{
					FillRect(topLeft_X, topLeft_Y, boxWidth, boxHeight, backColour);
					FillRect(topLeft_X + 1, topLeft_Y + 1, boxWidth - 2, boxHeight - 2, olc::GREEN);
				}
			}
		}

		//solveDijkstra();
		if (GetKey(olc::Key::K1).bPressed)
		{
			solveAStar();
			bRefreshPath = true;
		}
		if (GetKey(olc::Key::K2).bPressed)
		{
			solveDijkstra(); 
			bRefreshPath = true;
		}

		if (GetKey(olc::Key::K0).bPressed) 
		{

			bool bFinished = false;
			while (!bFinished) 
			{
				if (!recursiveStack.empty()) {
					recursiveStack.pop();
				}
				else 
				{
					bFinished = true;
				}
			}
			
			bRefreshPath = false;

			srand(clock());
			createMaze();
		}

		if (bRefreshPath) {
			if (nodeEnd != nullptr)
			{
				sNode* p = nodeEnd;
				while (p->parent != nullptr)
				{
					DrawLine(p->x * boxWidth + boxWidth / 2, p->y * boxHeight + boxHeight / 2, p->parent->x * boxWidth + boxWidth / 2, p->parent->y * boxHeight + boxHeight / 2, olc::YELLOW);
					// Set next node to this node's parent
					p = p->parent;
				}
			}
		}
		return true;
	}
};

int main()
{
	//seed for random number generator
	srand(clock());

	pathFindingAlgorithm demo;
	if (demo.Construct(SCREENHEIGHT, SCREENWIDTH, 1, 1))
		demo.Start();
	return 0;  
}