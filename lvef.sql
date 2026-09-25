-- This code file extracts LVEF from the mimiciv_echo.echo table. Note: the extracted lvef is text, not numeric, and must be converted to numeric after importing into R.

-- First extract all records related to LVEF
DROP TABLE IF EXISTS public.lvef1; 
CREATE TABLE public.lvef1 AS
SELECT *
FROM mimiciv_echo.echo
WHERE measurement IN ('lvef_3d','biplane_lvef','rest_biplane_lvef','lvef','rest_lvef')
AND result IS NOT NULL;
SELECT COUNT(*) FROM public.lvef1; -- 192351 records

-- For the same patient (subject_id) at the same examination time (measurement_datetime), when multiple LVEF measurement methods coexist, keep only one, according to the priority of the above options from left to right.
DROP TABLE IF EXISTS public.lvef; 
CREATE TABLE public.lvef AS
SELECT DISTINCT ON (subject_id, measurement_datetime)
    subject_id,
    measurement_id,
    measurement_datetime,
    measurement,
    result AS lvef,
    unit
FROM public.lvef1
ORDER BY 
    subject_id,
    measurement_datetime,
    array_position(ARRAY['lvef_3d','biplane_lvef','rest_biplane_lvef','lvef','rest_lvef'], measurement);
SELECT COUNT(*) FROM public.lvef; -- 191351 records
DROP TABLE IF EXISTS public.lvef1; 
