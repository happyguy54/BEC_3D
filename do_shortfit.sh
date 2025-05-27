#!/bin/bash
# Function to run a command and handle errors
run_command() {
  echo "Running: $@"
  "$@"
  if [ $? -ne 0 ]; then
    echo "ERROR: Command failed: $@"
    exit 1
  fi
}


rej_osl_args=""
indices=""
fit_option=""
fit_value=""
pT=""
qMin=0.02
QMax=1.82
r2=false

# Parse arguments dynamically
while [[ $# -gt 0 ]]; do
  case "$1" in
    --r2)
      r2=true
      shift
      ;;
    -c2)
      indices=$2
      shift 2
      ;;
    -ref)
      fit_option=$1
      fit_value=$2
      shift 2
      ;;
    -pt)
      pT=$2
      shift 2
      ;;
    -qmin)
      qMin=$2
      shift 2
      ;;
    -qmax)
      QMax=$2
      shift 2
      ;;
    --rej-fin)
      # Append --rej-fin and the contents of rej_OSL_ULS.txt
      rej_osl_args+="--rej-fin "
      shift
      # end the loop
      ;;
    *)
      rej_osl_args+=$@
      break
      ;;
  esac
done

# Check if the fit_option argument is "-ref" and fit_value is "ohp"
if [ "$fit_option" == "-ref" ] && [ "$fit_value" == "ohp" ]; then
  fit_cmd="-o"
elif [ "$fit_option" == "-ref" ] && [ "$fit_value" == "mix" ]; then
  fit_cmd="-m"
else
  fit_cmd="-u"
fi

for col in "alpha_out" "alpha_side" "alpha_long"; do
  declare $col=""
done
alpha_cmd=""
echo "Indices: $indices"  # Debug line
echo "rej_OSL arguments: $rej_osl_args"  # Debug line
echo "Fit command: $fit_cmd"  # Debug line

# Split the indices string into an array
IFS=',' read -ra IDX <<< "$indices"
IFS=',' read -ra qMaxval <<< "$QMax"
IFS=',' read -ra pt <<< "$pT"

echo "IDX: ${IDX[@]}"  # Debug line
echo "qMaxval: ${qMaxval[@]}"  # Debug line
echo "pT: ${pt[@]}"  # Debug line

# Read the YAML file
yaml_file="input/input.yaml"
temp_yaml_file="input/temp_input.yaml"

# Initialize variables
noncom=true
inside_uls_section=false
inside_ohp_section=false
inside_mix_section=false
lines=0
file_list=()
# Read the YAML file line by line
while IFS= read -r line || [ -n "$line" ]; do
  # Check if we are in the "mix" section
  if [ "$fit_cmd" == "-u" ]; then
    if [[ $line =~ ^uls: ]]; then
      inside_uls_section=true
      continue
    fi
  fi
  if [ "$fit_cmd" == "-o" ]; then
    if [[ $line =~ ^ohp: ]]; then
      inside_ohp_section=true
      continue
    fi
  fi
  if [ "$fit_cmd" == "-m" ]; then
    if [[ $line =~ ^mix: ]]; then
      inside_mix_section=true
      continue
    fi
  fi

  if $inside_uls_section; then
    if [ "$lines" -lt 2 ]; then
      # Extract the file path and add it to the list (including commented lines)
      file_path=$(echo "$line" | sed -E "s/.*'(\/media\/[^']+)'.*/\1/")
      file_list+=("$file_path")
      lines=$((lines + 1))
      # echo "Found file path: $file_path"  # Debug line
    fi
  fi
  if $inside_ohp_section; then
    if [ "$lines" -lt 4 ]; then
      # Extract the file path and add it to the list (including commented lines)
      file_path=$(echo "$line" | sed -E "s/.*'(\/media\/[^']+)'.*/\1/")
      file_list+=("$file_path")
      lines=$((lines + 1))
      # echo "Found file path: $file_path"  # Debug line
    fi
  fi
  # If inside the "mix" section, extract file paths starting with "/media/"
  if $inside_mix_section; then
    if [[ $line =~ /media/ ]]; then
        if $noncom && [[ $line =~ ^[[:space:]]*# ]]; then
          continue
        fi
      # Extract the file path and add it to the list (including commented lines)
      file_path=$(echo "$line" | sed -E "s/.*'(\/media\/[^']+)'.*/\1/")
      file_list+=("$file_path")
      # echo "Found file path: $file_path"  # Debug line
    fi
  fi
done < "$yaml_file"
# echo "File list: ${file_list[@]}"  # Debug line
# Loop through the pt values and create a temporary YAML content with the replaced pt value
if [ -z "$pT" ]; then
  pt=("100")
fi

# for pT_value in "${pt[@]}"; do
#   if [ "$pT_value" == "default" ]; then
#     cp "$yaml_file" "$temp_yaml_file"
#   else
#     sed "s/pt_[^_]*_/pt_${pT_value}_/g" "$yaml_file" > "$temp_yaml_file"
#   fi
for pT_value in "${pt[@]}"; do
  # Process the file list in pairs if 
  
  for ((k = 0; k < ${#file_list[@]}; k += 2)); do
    # Get the first and second file in the pair
    if $r2; then
      file1="${file_list[k]}"
      file2="${file_list[k+1]}"
      file3="${file_list[k+2]}"
      file4="${file_list[k+3]}"

      file1_replaced=$(echo "$file1" | sed "s/pt_[0-9]\+_trk/pt_${pT_value}_trk/")
      file2_replaced=$(echo "$file2" | sed "s/pt_[0-9]\+_trk/pt_${pT_value}_trk/")
      file3_replaced=$(echo "$file3" | sed "s/pt_[0-9]\+_trk/pt_${pT_value}_trk/")
      file4_replaced=$(echo "$file4" | sed "s/pt_[0-9]\+_trk/pt_${pT_value}_trk/")

      cat > "$temp_yaml_file" <<EOL
mix:
  file_data: '$file1_replaced'
  file2_data: '$file2_replaced'
  file_mc: '$file3_replaced'
  file2_mc: '$file4_replaced'
EOL
    echo "Generated temp YAML for r2=true with pt_value $pT_value:"
    cat "$temp_yaml_file"
    k=$((k+2))
    else
      file1="${file_list[k]}"
      file2="${file_list[k+1]}"

      # Replace the "pt_x" placeholder with "pt_${pT_value}" in the file paths
      file1_replaced=$(echo "$file1" | sed "s/pt_[0-9]\+_trk/pt_${pT_value}_trk/")
      file2_replaced=$(echo "$file2" | sed "s/pt_[0-9]\+_trk/pt_${pT_value}_trk/")

      # Determine if the pair is data or MC based on the file path
      if [[ $file1 =~ /data/ ]]; then
        pair_type="data"
        cat > "$temp_yaml_file" <<EOL
mix:
  file_data: '$file1_replaced'
  file2_data: '$file2_replaced'
EOL
      elif [[ $file1 =~ /mc/ ]]; then
        pair_type="mc"
        cat > "$temp_yaml_file" <<EOL
mix:
  file_mc: '$file1_replaced'
  file2_mc: '$file2_replaced'
EOL
      fi
    fi

    # Display the generated temp YAML file
    echo "Generated temp YAML for $pair_type with pt_value $pT_value:"
    cat "$temp_yaml_file"

    # Get today's date in the format used in the folder structure
    today=$(LC_TIME=C date +%Y_%b_%d)

    # Find the latest created folder (c2_data, c2_mc, or r2) in ./output/today's_date/
    

    echo "Using prefix: $prefix"
    # Loop over the indices and call the commands
    for i in "${IDX[@]}"; do
      for qMax in "${qMaxval[@]}"; do
        echo "Running commands with index $i and pt_${pT_value}"
        if [ "$i" -eq 4 ]; then
          latest_csv=""
        else
          # latest_csv=$(ls -t ./output/"$today"/"$prefix"/csv/*mult_proj*.csv | head -n 1)
          latest_folder=$(ls -td ./output/"$today"/*/ | head -n 1)

          # Extract the folder name (e.g., c2_data, c2_mc, or r2)
          if [ -n "$latest_folder" ]; then
            prefix=$(basename "$latest_folder")
          else
            echo "Error: No folder found in ./output/$today/"
            exit 1
          fi

          latest_csv=$(ls -t ./output/"$today"/"$prefix"/csv/*mult_proj*.csv | head -n 1) #in case of c2 being 4,6,8,10
          #   latest_csv=$(ls -t ./output/*/c2_${i}/csv/*mult_proj*.csv | head -n 1)
          # else
          #   latest_csv=$(ls -t ./output/*/c2_$((i-2))/csv/*mult_proj*.csv | head -n 1) #in case of c2 being 4,6,8,10
          # fi
          if [ -z "$latest_csv" ]; then
            echo "Error: No CSV file found"
            exit 1
          else
            echo "Latest CSV file: $latest_csv"
          fi
          # Get the column number of #alpha_{out} from the first line
          for col in "#alpha_{out}" "#alpha_{side}" "#alpha_{long}"; do
            col_name=$(echo $col | tr -d '#{}')
            col_num=$(awk -F, -v col="$col" 'NR==1{for(i=1;i<=NF;i++){if($i==col){print i;exit}}}' $latest_csv)
            value=$(LC_NUMERIC=C awk -v col=$col_num -F, '$1 == 21 && $2 == 250 {printf "%.2f", $col}' $latest_csv)
            echo "Value of $col_name: $value"
            eval "$col_name=$value"
          done
        fi
        # alpha_out=1.16
        # alpha_side=2.5
        # alpha_long=1.29
        alpha_cmd=""
        #make if for ohp and mix
        if [ "$fit_cmd" == "-o" ] || [ "$fit_cmd" == "-m" ]; then
          if [ "$i" -gt 4 ]; then
            alpha_cmd+=" --alpha-out "$alpha_out
          fi
          if [ "$i" -gt 6 ]; then
            alpha_cmd+=" --alpha-long "$alpha_long
          fi
          if [ "$i" -gt 8 ]; then
            alpha_cmd+=" --alpha-side "$alpha_side
          fi
        fi
        if [ "$fit_cmd" == "-u" ]; then
          if [ "$i" -gt 4 ]; then
            alpha_cmd+=" --alpha-long "$alpha_long
          fi
          if [ "$i" -gt 6 ]; then
            alpha_cmd+=" --alpha-side "$alpha_side
          fi
          if [ "$i" -gt 8 ]; then
            alpha_cmd+=" --alpha-out "$alpha_out
          fi
        fi
        #sourceroot
        echo "./do_fit.py $fit_cmd $alpha_cmd --unf -n --input-yaml input/temp_input.yaml --c2-index $i $rej_osl_args --q-min $qMin --q-max $qMax"
        run_command ./do_fit.py $fit_cmd $alpha_cmd --unf -n --input-yaml input/temp_input.yaml --c2-index $i $rej_osl_args --q-min $qMin --q-max $qMax
        echo "./do_fit.py $fit_cmd $alpha_cmd -k --input-yaml input/temp_input.yaml --c2-index $i $rej_osl_args --q-min $qMin --q-max $qMax"
        run_command ./do_fit.py $fit_cmd $alpha_cmd -k --input-yaml input/temp_input.yaml --c2-index $i $rej_osl_args --q-min $qMin --q-max $qMax
        echo "./do_fit.py $fit_cmd $alpha_cmd --tt --input-yaml input/temp_input.yaml --c2-index $i $rej_osl_args --q-min $qMin --q-max $qMax"
        run_command ./do_fit.py $fit_cmd $alpha_cmd --tt --input-yaml input/temp_input.yaml --c2-index $i $rej_osl_args --q-min $qMin --q-max $qMax
        echo "./do_fit.py $fit_cmd $alpha_cmd --input-yaml input/temp_input.yaml --c2-index $i $rej_osl_args --q-min $qMin --q-max $qMax"
        run_command ./do_fit.py $fit_cmd $alpha_cmd --input-yaml input/temp_input.yaml --c2-index $i $rej_osl_args --q-min $qMin --q-max $qMax
      done
    done
  done
done