#include "../include/main.h"

/*
An implementation of the Fiduccia-Mattheyses partitioning algorithm
Jan Kowalski 3/2020-8/2020
*/

// Main control function for the program
// The ibm_testbench_number is 1-18, 0 for files not in the testbench suite.
void import_data_and_run_algorithm(char *are_filename, char *netD_filename)
{
	clock_t begin = clock();

	// Obtain information from the two data files (are, netD)
	struct condensed *information = read_in_data_to_arrays(are_filename, netD_filename);
	// Add useful information about partition sizes
	information->desired_area = (int)(RATIO * information->total_area);
	information->ratio = RATIO;

	// Set the FM_chromosome to NULL so that GA proceeds normally
	information->FM_chromosome = NULL;

	int i;
	int lowest_global_cutsize = 99999999;
	long area_util = 0;	 // 存每輪最佳解
	long area2_util = 0; // 存每輪最佳解
	for (i = 0; i < FM_NUM_PASSES; i++)
	{

		if (i > 0)
		{
			reset_cells_and_nets(information);
			delete_partition(information->partition_A);
			delete_partition(information->partition_B);
		}
		// Create and initialize the two partition structs
		initialize_two_partitions(information);

		// Separate the cells into one of the two partitions
		populate_partitions(information);
		////////////////////print cell output
		/*
		FILE *Output = fopen("Output_initial", "w");
		fprintf(Output, "Lowest cutstate achieved: %d\n", information->lowest_cutstate);
		fprintf(Output, "A : %ld , B : %ld\n", information->initial_area, information->initial_area2);
		fprintf(Output, "max_in_A : %d , max_in_B : %d\n", information->tolerance, information->tolerance1);
		double uA = (double)information->die_area * information->utilA / 100;
		double uB = (double)information->die_area * information->utilA / 100;
		fprintf(Output, "die_size : %ld maxutil_area_A : %f , maxutil_area_B : %f\n", information->die_area, uA, uB);
		int area0;
		int area1;
		struct node *node0 = information->partition_A->cells_in_partition->head->next;
		while (node0 != information->partition_A->cells_in_partition->tail)
		{
			struct cell *cell0 = (struct cell *)node0->data_structure;
			area0 = 0;
			if (cell0->which_partition == 0)
				area0 = cell0->area;
			fprintf(Output, "%d      %d           %4d    %4d\n", cell0->identifier + 1, cell0->which_partition, cell0->area, cell0->area2); //
			node0 = node0->next;
		}
		struct node *node1 = information->partition_B->cells_in_partition->head->next;
		while (node1 != information->partition_B->cells_in_partition->tail)
		{
			struct cell *cell1 = (struct cell *)node1->data_structure;
			area1 = 0;
			if (cell1->which_partition == 1)
				area1 = cell1->area2;
			fprintf(Output, "%d      %d           %4d    %4d\n", cell1->identifier + 1, cell1->which_partition, cell1->area, cell1->area2); //
			node1 = node1->next;
		}
		fclose(Output);
		*/
		////////////////////
		// Reset the partition information with new GA info if FM begins to climb
		if ((information->FM_chromosome != NULL) && (information->FM_chromosome->cutstate > (lowest_global_cutsize + GA_TRIGGER)))
		{
			delete_partition(information->partition_A);
			delete_partition(information->partition_B);
			initialize_two_partitions(information);
			segregate_cells_with_GA(information);
		}

		// Important to create FM_chromosome after populate_partitions, as it checks for NULL to determine if GA or not
		// FM_chromosome is no longer needed, as its information has been used to population partitions

		delete_chromosome(information->FM_chromosome);

		information->FM_chromosome = malloc(sizeof(struct chromosome));
		initialize_chromosome(information->FM_chromosome, information);
		printf("%d, %d\n", information->tolerance, information->tolerance1);
		// Run the algorithm
		// printf("%d, %d\n", information->partition_A->total_partition_area, information->partition_B->total_partition_area);
		fiduccia_mattheyses_algorithm(information);

		if (information->lowest_cutstate < lowest_global_cutsize)
		{
			lowest_global_cutsize = information->lowest_cutstate;
			printf("lowest_global_cutsize : %d\n", lowest_global_cutsize);
			/////////存入每輪最佳解面積
			area_util = information->final_area;
			area2_util = information->final_area2;
			/////////
			// ////////////////////print cell output
			// FILE *Output = fopen("Output_true", "w");
			// printf("Output_true 開\n");
			// fprintf(Output, "Lowest cutstate achieved: %d, %d\n", information->lowest_cutstate, lowest_global_cutsize);
			// fprintf(Output, "A : %ld , B : %ld\n", area_util, area2_util);
			// fprintf(Output, "max_in_A : %d  Macro : %d , max_in_B : %d  Macro : %d\n", information->tolerance, information->tolerance_Macro, information->tolerance1, information->tolerance1_Macro);
			// double uA = (double)information->die_area * information->utilA / 100;
			// double uB = (double)information->die_area * information->utilA / 100;
			// fprintf(Output, "maxutil_A : %f , maxutil_B : %f\n", (float)information->final_area / information->die_area, (float)information->final_area2 / information->die_area);
			// fprintf(Output, "die_size : %ld maxutil_area_A : %f , maxutil_area_B : %f\n", information->die_area, uA, uB);
			// long area0 = 0;
			// long area1 = 0;
			// struct node *node0 = information->partition_A->cells_in_partition->head->next;
			// while (node0 != information->partition_A->cells_in_partition->tail)
			// {
			// 	struct cell *cell0 = (struct cell *)node0->data_structure;
			// 	if (cell0->which_partition == 0)
			// 		area0 += cell0->area;
			// 	if (cell0->which_partition == 1)
			// 		area1 += cell0->area2;
			// 	fprintf(Output, "%d      %d           %4d    %4d\n", cell0->identifier + 1, cell0->which_partition, cell0->area, cell0->area2); //
			// 	node0 = node0->next;
			// }
			// struct node *node1 = information->partition_B->cells_in_partition->head->next;
			// while (node1 != information->partition_B->cells_in_partition->tail)
			// {
			// 	struct cell *cell1 = (struct cell *)node1->data_structure;
			// 	if (cell1->which_partition == 1)
			// 		area1 += cell1->area2;
			// 	if (cell1->which_partition == 0)
			// 		area0 += cell1->area;
			// 	fprintf(Output, "%d      %d           %4d    %4d\n", cell1->identifier + 1, cell1->which_partition, cell1->area, cell1->area2); //
			// 	node1 = node1->next;
			// }
			// fprintf(Output, "%ld,%ld", area0, area1);
			// fclose(Output);
			// printf("Output_true 關\n");
			// ////////////////////
		}
		// Some information needs to be freed between repeats
		free(information->access_);
	}

	printf("FM_NUM_PASSES: %d\n", FM_NUM_PASSES);
	printf("Lowest cutstate achieved: %d, %d\n", information->lowest_cutstate, lowest_global_cutsize);
	float A = (float)/*(information->final_area)*/ area_util / information->die_area;
	float B = (float)/*(information->final_area2)*/ area2_util / information->die_area;
	printf("A : %ld,%ld , B : %ld,%ld\n", area_util, information->final_area, area2_util, information->final_area2);
	printf("maxutil_A : %f , maxutil_B : %f\n", A, B);
	////////////////////print cell output
	FILE *Output = fopen("Output", "w");
	fprintf(Output, "Lowest cutstate achieved: %d, %d\n", information->lowest_cutstate, lowest_global_cutsize);
	fprintf(Output, "A : %ld , B : %ld\n", area_util, area2_util);
	fprintf(Output, "max_in_A : %d  Macro : %d , max_in_B : %d  Macro : %d\n", information->tolerance, information->tolerance_Macro, information->tolerance1, information->tolerance1_Macro);
	double uA = (double)information->die_area * information->utilA / 100;
	double uB = (double)information->die_area * information->utilA / 100;
	fprintf(Output, "maxutil_A : %f , maxutil_B : %f\n", A, B);
	fprintf(Output, "die_size : %ld maxutil_area_A : %f , maxutil_area_B : %f\n", information->die_area, uA, uB);
	long area0 = 0;
	long area1 = 0;
	struct node *node0 = information->partition_A->cells_in_partition->head->next;
	while (node0 != information->partition_A->cells_in_partition->tail)
	{
		struct cell *cell0 = (struct cell *)node0->data_structure;
		if (cell0->which_partition == 0)
			area0 += cell0->area;
		if (cell0->which_partition == 1)
			area1 += cell0->area2;
		fprintf(Output, "%d      %d           %4d    %4d\n", cell0->identifier + 1, cell0->which_partition, cell0->area, cell0->area2); //
		node0 = node0->next;
	}
	struct node *node1 = information->partition_B->cells_in_partition->head->next;
	while (node1 != information->partition_B->cells_in_partition->tail)
	{
		struct cell *cell1 = (struct cell *)node1->data_structure;
		if (cell1->which_partition == 1)
			area1 += cell1->area2;
		if (cell1->which_partition == 0)
			area0 += cell1->area;
		fprintf(Output, "%d      %d           %4d    %4d\n", cell1->identifier + 1, cell1->which_partition, cell1->area, cell1->area2); //
		node1 = node1->next;
	}
	fprintf(Output, "%ld,%ld", area0, area1);
	fclose(Output);
	////////////////////
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	if (PRINT_EXECUTION_TIME)
		printf("Program execution time: %f\n", time_spent);

	free_all_memory(information);
}

// File management, options for demo
int main()
{

	// Run the demo if the user option is enabled
	if (RUN_DEMO_WITH_TESTDATA)
	{
		printf("################################################\n");
		printf("Running demo with visualized partitions\n");
		printf("Cells are visible only if they can switch sides\n");
		printf("################################################\n");

		import_data_and_run_algorithm(TEST_ARE_FILENAME, TEST_NETD_FILENAME);

		printf("######################################################\n");
		printf("Conclusion of demo. To disable, access include/main.h\n");
		printf("######################################################\n");

		return 0;
	}

	// Define format for ibm benchmark files
	char are_filename[22] = "data/ibm00.are";
	char netD_filename[22] = "data/ibm00.netD";

	char number[3];
	snprintf(number, 10, "%d", IBM_FILE_NUMBER);
	if (IBM_FILE_NUMBER < 10)
	{
		memcpy(are_filename + 9, number, 1);
		memcpy(netD_filename + 9, number, 1);
	}
	else
	{
		memcpy(are_filename + 8, number, 2);
		memcpy(netD_filename + 8, number, 2);
	}

	struct stat buffer;
	FILE *file_check1, *file_check2;
	// Check if files are available
	if ((stat(are_filename, &buffer) == 0) && (stat(netD_filename, &buffer) == 0))
	{
		import_data_and_run_algorithm(are_filename, netD_filename);
	}
	else
	{
		printf("Either %s or %s is inaccessible by the program\n", are_filename, netD_filename);
	}

	return 0;
}

void reset_cells_and_nets(struct condensed *information)
{
	int i;
	struct cell *temp_cell;
	struct net *temp_net;
	struct node *locked_start;
	struct node *locked_end;
	struct node *free_end;

	for (i = 0; i < information->CELL_array_size; i++)
	{
		temp_cell = information->CELL_array[i];
		temp_cell->cell_state = FREE;
		temp_cell->gain = 0;
		temp_cell->GAIN_array_node = NULL;
	}

	for (i = 0; i < information->NET_array_size; i++)
	{
		// O(1) transfer of all locked cells to each free_cells list
		temp_net = information->NET_array[i];
		locked_start = temp_net->locked_cells->head->next;
		free_end = temp_net->free_cells->tail->previous;
		connect_two_nodes(free_end, locked_start);
		locked_end = temp_net->locked_cells->tail->previous;
		connect_two_nodes(locked_end, temp_net->free_cells->tail);
		connect_two_nodes(temp_net->locked_cells->head, temp_net->locked_cells->tail);

		temp_net->num_cells_in_[PARTITION_A] = 0;
		temp_net->num_cells_in_[PARTITION_B] = 0;
	}

	information->lowest_cutstate = information->current_cutstate;
}

void free_all_memory(struct condensed *information)
{
	int i;
	// Partitions should be deallocated before nets
	delete_partition(information->partition_A);
	delete_partition(information->partition_B);
	// Nets should be deallocated before cells
	for (i = 0; i < information->NET_array_size; i++)
	{
		delete_net(information->NET_array[i]);
	}
	// Cells deallocated last
	for (i = 0; i < information->CELL_array_size; i++)
	{
		delete_cell(information->CELL_array[i]);
	}
	delete_chromosome(information->FM_chromosome);
	free(information->NET_array);
	free(information->CELL_array);
	free(information);
}
