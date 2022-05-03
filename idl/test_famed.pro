pro test_famed

catalog_id = 'KIC'
star_id = '012069424'
teff = 5825

; Do some preliminary checks on the directories for Background (in case the related software is not installed)
if file_test('../../Background/data',/directory) eq 0 then begin
	file_mkdir,'../../Background/data'
endif	

if file_test('../../Background/results',/directory) eq 0 then begin
	file_mkdir,'../../Background/results'
endif	

; Remove existing test files if present
if file_test('../../PeakBagging/results/' + catalog_id + star_id + '/*') eq 1 then begin
	print,'Do you want to remove old test files? ( 1 = YES | Otherwise = NO)'
	read,response
	if response eq 1 then begin
		if file_test('../../PeakBagging/results/' + catalog_id + star_id + '/old',/directory) eq 0 then begin
			file_mkdir,'../../PeakBagging/results/' + catalog_id + star_id + '/old'
		endif else begin
			spawn,'rm -r ' + '../../PeakBagging/results/' + catalog_id + star_id + '/old/*',output,/stderr
		endelse
		
		spawn,'mv ../../PeakBagging/results/' + catalog_id + star_id + '/* ../../PeakBagging/results/' + catalog_id + star_id + '/old/',output,/stderr
	endif else begin
		print,'Previous test files were not deleted. Aborting test.'
		return
	endelse
endif

; Copy stellar PSD and background fit solution inside the expected directories
spawn,'cp ../tutorials/data/Background/data/' + catalog_id + star_id + '.txt ../../Background/data/',output,/stderr 
spawn,'cp -r ../tutorials/data/Background/results/' + catalog_id + star_id + ' ../../Background/results/',output,/stderr
start_famed,catalog_id,star_id,teff,/fit,/global,/chunk
end

pro test_external_background

catalog_id = 'KIC'
star_id = '012008916'
teff = 5454

spawn,'cp famed_configuring_parameters.txt famed_configuring_parameters.txt.local'
spawn,'sed -i.old "22s^-99^../tutorials/data/Background/results/' + catalog_id + star_id + '/^g" famed_configuring_parameters.txt'

; Remove existing test files if present
if file_test('../../PeakBagging/results/' + catalog_id + star_id + '/*') eq 1 then begin
	print,'Do you want to remove old test files? ( 1 = YES | Otherwise = NO)'
	read,response
	if response eq 1 then begin
		if file_test('../../PeakBagging/results/' + catalog_id + star_id + '/old',/directory) eq 0 then begin
			file_mkdir,'../../PeakBagging/results/' + catalog_id + star_id + '/old'
		endif else begin
			spawn,'rm -r ' + '../../PeakBagging/results/' + catalog_id + star_id + '/old/*',output,/stderr
		endelse

		spawn,'mv ../../PeakBagging/results/' + catalog_id + star_id + '/* ../../PeakBagging/results/' + catalog_id + star_id + '/old/',output,/stderr
	endif else begin
		print,'Previous test files were not deleted. Aborting test.'
		spawn,'mv famed_configuring_parameters.txt.local famed_configuring_parameters.txt'
		return
	endelse
endif

start_famed,catalog_id,star_id,teff,/fit,/global,/chunk
spawn,'mv famed_configuring_parameters.txt.local famed_configuring_parameters.txt'

end