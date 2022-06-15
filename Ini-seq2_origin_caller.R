##G.Guilbaud JE Sale group 2021
##Ori peak caller and effiency score script

#This script uses R command line to (mainly) run bedtools and awk commands
#After counting the number of sequencing reads by windows of 100bp  (size can be defined by user), the number of reads per window is normalised to the total number of reads for each reaction.





###  Required software and Input file

#Input: aligned read files in bed format labelled with prefix 'Ini_HL_' or 'Ini_LL' and placed in the Input folder. For one Ini-seq experiment, the same suffix should be used for HL and LL.
#Segmented genome of the desired window: Files to be placed in './Ref.genome folder'. Window size must be stated at the end of filename eg: 'hg38_chr_size_windows_of_100.bed'
#Required to run this script: a subfolder './Output', itself containing two folders named as: './HL' and './LL'.




### Settings

			##Windows					

	Windows <-  c(100) ##several window sizes can be specified. They will be sequentially processed, one segmented genome file per desired window should be present in the subfolder Ref.genome.	
	
			##Minimum (Normalised) Tags considered for the analysis:
				#Factor driving threshold under which number of tags are changed by 'Value.HL.is.changed.by' & 'Value.LL.is.changed.by'. 
				#Set as a function of total number of reads (eg 50e6 reads, want to change ,5tags=0 => tags_factor=10.e-8. Previously HL_min_tags_factor <- 15.e-8, LL_min_tags_factor <- 5e-8)
	
	HL_min_tags_factor <- 13.e-8 
	LL_min_tags_factor <- 5.e-8
	
	Value.HL.is.changed.by <- 0
	Value.LL.is.changed.by <- 1e-20 ##LL value is changed by 1e-20 if equal to 0 or if lower than LL_min_tags_factor (thus avoiding dividing by 0)


	
			##Min log2 threshold for a window to be retained after counting reads per window and log2(HL/LL) computed

	log.threshold <- 2
	
	
			##Size (in bp) of adjacent retained windows to merge (the merged windows being later called an island.)

	Merge.adjacent.island <- 510
	
	
			## Final minimum island size

	Final.min.lim.Island.size <- 200	
	
	
			##Final enrichment threshold. (Once islands determined, recount reads in each island, compute HL/(HL+LL), and keep only island with ratio > HL.enrichment.cut.off)

	HL.enrichment.cut.off <- 0.8
	


### Remove temp files from previous run (if not yet done) 

		##Delete temp files and previous versions 
		
delete.temp.files <- paste("rm Input/*_single_nt_position.bed")
		try(system(delete.temp.files))	
		
delete.Output.files <- paste("rm Output/*.bed")
		try(system(delete.Output.files))
			
delete.Output.files.HL <- paste("rm Output/HL/*.bed")
		try(system(delete.Output.files.HL))	
			
delete.Output.files.LL <- paste("rm Output/LL/*.bed")
		try(system(delete.Output.files.LL))	
		
	
		
		
### Load and check Prefix of Input files  
			
File.list.Input <- list.files("Input/")[grep("*.bed",list.files("Input/"))]   ##list bed files only
	
File.prefix <- substr(File.list.Input,0,6)	##Keep first 6 characters of file name
		
		##Check PRefix of input file, stop script if not correct
		
for(prefix in File.prefix)	{
	
	if(prefix!="Ini_HL"){ if (prefix!="Ini_LL"){stop("Files prefix Have to be labelled as: Ini_HL_XXX.bed or Ini_LL_XXX_.bed")}}
	
}	
		
		
		
	
###  Loop for file and window(s)     	  			  

Unique.File.name <- unique(substr(File.list.Input,8,nchar(File.list.Input))) ##Retreive file name w/o Ini_XL (Allows several Ini-seq conditions to be run at once)

	for (Input.name in Unique.File.name ){ ##Loop to run experiemnts 
		
		for (Windows.size in Windows){	##Work on each window sequntially
			
			for (Conditon.to.run in c(1,2)){
			
			if(Conditon.to.run==1){Input.bed <- paste("Ini_HL_", Input.name,sep="")} #Re-asign HL prefix 
			if(Conditon.to.run==2){Input.bed <- paste("Ini_LL_", Input.name,sep="")} #Re-asign LL prefix 

		
		print(paste("Working on: ", Input.bed,sep=""))
		

			
	File.list <- list.files("Ref.genome/")[grep("*.txt",list.files("Ref.genome/"))]   ##list ref genome (only .txt), must match windows size selected


	
###  Compute middle position of each read and pass it to the start/end of each read


			try(system(paste("wc -l Input/",Input.bed," > Output/Number.reads.in.", Input.bed ,sep=""))) ###Count number of lines == number of reads
		

			Total.num.mapped.read <- read.table(paste("Output/Number.reads.in.", Input.bed ,sep=""))[1] ### Retreive number of lines 
		
			
			File.Prefix <- substr(Input.bed,0,nchar(Input.bed)-4) ##rm ".bed"

	
			##Call system to run awk command => Set up at 25 (for 50bp reads).			
	     	try(system( paste("awk '{print $1, $2+35, $2+35}' OFS='\t' Input/", Input.bed, " > ", "./Input/",File.Prefix,"_single_nt_position.bed",sep=""))) ##asign middle position. To be changed if reads length =! 50



	
### Count number of tags / windows


		
		print(paste("Counting tags in: ", Input.bed,sep=""))	
			

		Name.file.to.count <- paste(File.Prefix,"_single_nt_position.bed",sep="")

					
 			if(Conditon.to.run==1){	##run HL

 				Output.file.name <- paste("Output/HL/",substr(Name.file.to.count,0,nchar(Name.file.to.count)-4),"_by_", Windows.size,"_bp_windows.bed",sep="") ##Prepare name for output file
 				 				
				command <- paste("bedtools window -c -w 0 -a",paste(" Ref.genome/",substr(File.list,0,nchar(File.list)-4),"_windows_of_", Windows.size,".bed",sep="")," -b Input/", Name.file.to.count," > ", Output.file.name,	sep="") ##Prepare bedtools command to count reads
		
				try(system(command))		##Call bedtool

					### Normalised all conditions by total number of reads 
																		
 					print(paste("Normalising: ", Name.file.to.count,", by: ",Total.num.mapped.read,sep=""))

			print(paste("Changing #tags in HL by 0 if less than ", round(Total.num.mapped.read*HL_min_tags_factor,0) ," tags",sep=""))
			
			###For HL: If number of tags < @ min (defined by user), replace by 0
				Temp.norm.HL <- read.table(Output.file.name)
				num <- which(Temp.norm.HL[,4]<as.numeric(round(Total.num.mapped.read*HL_min_tags_factor,0)))
				Temp.norm.HL[num,4] <- Value.HL.is.changed.by

				
				write.table(Temp.norm.HL, Output.file.name,col.names=F,row.names=F,quote=F)

				awk.command <- paste("awk '{print $1, $2, $3, $4/",Total.num.mapped.read,"}' ", Output.file.name," > ","Output/HL/Normalised_",substr(Name.file.to.count,0,nchar(Name.file.to.count)-4),"_by_", Windows.size,"_bp_windows.bed" , sep="") #Normalised by total number of reads

				try(system(awk.command))	##Call awk

		
		}
			

			
			if(Conditon.to.run==2){	##run LL

 				Output.file.name <- paste("Output/LL/",substr(Name.file.to.count,0,nchar(Name.file.to.count)-4),"_by_", Windows.size,"_bp_windows.bed",sep="") ##Prepare name for output file
 				 				
				command <- paste("bedtools window -c -w 0 -a",paste(" Ref.genome/",substr(File.list,0,nchar(File.list)-4),"_windows_of_", Windows.size,".bed",sep="")," -b Input/", Name.file.to.count," > ", Output.file.name,	sep="") ##Prepare bedtools command
		
				try(system(command))	  ##Call bedtool


	##For LL: If number of tags < @ min (defined by user), replace by 0
	
				Temp.norm.LL <- read.table(Output.file.name)
				num <- which(Temp.norm.LL[,4]==0 & Temp.norm.LL[,4]<= as.numeric(round(Total.num.mapped.read* LL_min_tags_factor,0)))
				Temp.norm.LL[num,4] <- Value.LL.is.changed.by
				
				write.table(Temp.norm.LL, Output.file.name,col.names=F,row.names=F,quote=F)

					###Normalised all conditions by total number of reads 
			
 									print(paste("Normalising: ", Name.file.to.count,", by: ",Total.num.mapped.read,sep=""))

		
				awk.command.zero <- paste("awk '{print $1, $2, $3, $4/",Total.num.mapped.read,"}' ", Output.file.name," > ","Output/LL/Normalised_",substr(Name.file.to.count,0,nchar(Name.file.to.count)-4),"_by_", Windows.size,"_bp_windows.bed" , sep="") #Normalised by total number of reads
				
				try(system(awk.command.zero))	 ##Call awk
		
		}	
} ##End Condition to run
								
									
									
														
	
### Do Log2 ratios & filter



Condition.unique.name <- paste(substr(Input.name,0,nchar(Input.name)-4),"_single_nt_position_by_", Windows.size,"_bp_windows.bed",sep="")
					
print(paste("Normalising HL/LL: ",Condition.unique.name,sep=""))

HL.file <- read.table(paste("Output/HL/Normalised_Ini_HL_",Condition.unique.name,sep=""),sep="")
LL.file <- read.table(paste("Output/LL/Normalised_Ini_LL_",Condition.unique.name,sep=""),sep="")
												
non.0 <- which(HL.file[,4]!=0) ##select windows with more than 0 tags 
Temp.tab.out <- HL.file
Temp.tab.out[non.0,4] <- round(log((HL.file[non.0,4]/LL.file[non.0,4]),2),3)  ##Do log ratio on HL not equal to 0
write.table(Temp.tab.out,paste("Output/For_IGV_Log_ratio_of_Ini_HL_LL",Condition.unique.name,sep=""),quote=F,col.names=F,row.names=F,sep="\t") ##write.table

	##Apply log.threshold
	Filtered.output <- Temp.tab.out[which(Temp.tab.out[,4]>= log.threshold),]

		##Need to remove Inf to please bedtool
		temp.pos <- which(Filtered.output[,4]>1000)
		Filtered.output[temp.pos,4] <- 10000  

write.table(Filtered.output,paste("Output/0.Filtered_Log_ratio_of_Ini_HL_LL_more_than_", log.threshold,"_of_",Condition.unique.name,sep=""),quote=F,col.names=F,row.names=F,sep="\t") 


### Merge adjacent windows  									 

	###Merge adjacent islands if < Merge.adjacent.island (i.e: adjacent ones!) 
		print(paste("Merging island ",sep=""))

command.2 <- paste("bedtools merge -d ", Merge.adjacent.island," -c 4 -o sum,count -i",paste(" Output/0.Filtered_Log_ratio_of_Ini_HL_LL_more_than_", log.threshold,"_of_",Condition.unique.name,sep="")," > Output/1.Merge_Filtered_Log_ratio_of_Ini_HL_LL_more_than_", log.threshold,"_of_",Condition.unique.name, sep="") ##Prepare bedtools command

		try(system(command.2))	##Call bedtool



### From called island, 're'-compute number of tags and normalised
		
		
		###Do awk on 1.Merge_Filtered_Log_ratio_of_Ini_HL_LL_more_than_", log.threshold,"_of_",Condition.unique.name => To keep only chr start end / as other info not useful 

	awk.command.clear.table <- paste("awk '{print $1, $2, $3}'  OFS='\t' ", paste("Output/1.Merge_Filtered_Log_ratio_of_Ini_HL_LL_more_than_", log.threshold,"_of_",Condition.unique.name, sep="") ," > ","Output/Tidy.of.table.1.temp.file.bed" , sep="")

				try(system(awk.command.clear.table))	##Call awk


		##Now open HL and LL *.gDNA_rm_dup_single_nt_position 
	
 File.list.Input.sngl.nt.pos <- list.files("Input/")[grep(paste(substr(Condition.unique.name,1,(nchar(Condition.unique.name)-nchar(paste("_by_",Windows.size,"_bp_windows.bed",sep="")))),sep=""),list.files("Input/"))]   ##list bed files only

		

for (Final.count in File.list.Input.sngl.nt.pos){
	
	
	###Split to count in eihter HL or LL
				HL <- 0
				LL <- 0
				
				HL <- length(grep("HL", Final.count))
				LL <- length(grep("LL", Final.count))
				
 			if(HL>0 & LL==0){	##run for HL and test if LL really = 0 


		##Prepare name for tab.out
		Output.file.name.final.count.HL <- paste("Output/HL/Recount.tags.per.called.island.in.HL.of.",Condition.unique.name,sep="")
	
		##Bedtool to count reads within selected regions
		command.3 <- paste("bedtools window -c -w 0 -a Output/Tidy.of.table.1.temp.file.bed -b Input/", Final.count," > ", Output.file.name.final.count.HL,	sep="") ##Prepare bedtools command to count reads
	
		print(paste("Recount_reads_in_called_region_of_",Final.count,sep=""))
		
		try(system(command.3))	##Call bedtools


	###Normalised number of counts by total number of reads 
			File.to.open <- list.files("Output")[grep(paste("Number.reads.in.",substr(Final.count,0,nchar(Final.count)-23),"*",sep=""),list.files("Output/"))]
			Total.num.mapped.read <- read.table(paste("Output/",File.to.open,sep=""))[1]
 									
 									print(paste("Normalising: ", Final.count,", by: ",Total.num.mapped.read,sep=""))

	##In this awk, col 4 = number of merged island and $5 number of tags, normalised col 5 by number of reads and swap 4 <-> 5
			awk.command.2 <- paste("awk '{print $1, $2, $3, $4/",Total.num.mapped.read,"}'  OFS='\t'  ", Output.file.name.final.count.HL," > ","Output/HL/Final.Normalised_Recount.tags.per.called.island.in.HL.of.",Condition.unique.name,sep="")
				
				try(system(awk.command.2))	##Call awk


	}
	
	
	
	if(LL>0 & HL==0){	##run for LL and test if HL really = 0 

	##Prepare name for tab.out
	
	Output.file.name.final.count.LL <- paste("Output/LL/Recount.tags.per.called.island.in.LL.of.",Condition.unique.name,sep="")
	
	##Bedtool to count reads within selected regions
	
	command.3 <- paste("bedtools window -c -w 0 -a Output/Tidy.of.table.1.temp.file.bed -b Input/", Final.count," > ", Output.file.name.final.count.LL,	sep="") ##Prepare bedtools command to count reads
	
										print(paste("Recount_reads_in_called_region_of_",Final.count,sep=""))

			try(system(command.3))	##Call bedtool
			
			

###Normalised number of counts by total number of reads 


			File.to.open <- list.files("Output")[grep(paste("Number.reads.in.",substr(Final.count,0,nchar(Final.count)-23),"*",sep=""),list.files("Output/"))]
			Total.num.mapped.read <- read.table(paste("Output/",File.to.open,sep=""))[1]
 									
 									print(paste("Normalising: ", Final.count,", by: ",Total.num.mapped.read,sep=""))

				##In this awk, col 4 = number of merged island and $5 number of tags, normalised col 5 by number of reads and swap 4 <-> 5

				awk.command.2 <- paste("awk '{print $1, $2, $3, $4/",Total.num.mapped.read,"}'  OFS='\t' ", Output.file.name.final.count.LL," > ","Output/LL/Final.Normalised_Recount.tags.per.called.island.in.LL.of.",Condition.unique.name,sep="")
				
				try(system(awk.command.2))	##Call awk

	}
}

delete.temp.files <- paste("rm Output/Tidy.of.table.1.temp.file.bed")
		try(system(delete.temp.files))	
		

### Do HL/(HL+LL) & filter 

					print(paste("Do HL/(HL+LL): ",Condition.unique.name,sep=""))

	HL.file <- read.table(paste("Output/HL/Final.Normalised_Recount.tags.per.called.island.in.HL.of.",Condition.unique.name,sep=""),sep="")
	LL.file <- read.table(paste("Output/LL/Final.Normalised_Recount.tags.per.called.island.in.LL.of.",Condition.unique.name,sep=""),sep="")

													
	##This should not be required as already selected region enriched in HL - Do test

	non.0 <- which(HL.file[,4]!=0)
	if(length(which(HL.file[,4]==0))>0){print("WARNING: REGIONS WITH 0 TAGS IN HL")}

		HL.file[,4] <- round(HL.file[,4]/(HL.file[,4]+LL.file[,4]),4)  

###Select only island > Final.min.lim.Island.size

Filtered.output.temp <- HL.file[which((HL.file[,3]-HL.file[,2])>= Final.min.lim.Island.size),]


###Threshold HL/(HL+LL) (HL.enrichment.cut.off)

Filtered.output <- Filtered.output.temp[which(Filtered.output.temp[,4]>=HL.enrichment.cut.off),]
	
		
write.table(Filtered.output,paste("Output/Final_HL_enrchiment_cutoff_", HL.enrichment.cut.off,"log_threshold_called_island_of_", log.threshold,"_for_", Condition.unique.name,sep=""),quote=F,col.names=F,row.names=F,sep="\t") ##write.table





	}##End of windows loop
	
} ##End Input.name loop

	print("Done")
