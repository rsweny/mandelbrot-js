function matrix(rows, cols, defaultValue) {
	var arr = [];
	for(var i = 0; i < rows; i++) {
		arr.push([]);
		arr[i].push(new Array(cols));

		for(var j=0; j < cols; j++) {
			arr[i][j] = defaultValue;
		}
	}
	return arr;
}

function RotateX(angle)
{
	var rot = matrix(3,3,0.0);
	var s = Math.sin(angle);
	var c = Math.cos(angle);
	rot[0][0] = 1.0;
	rot[1][1] = c;
	rot[2][2] = c;
	rot[1][2] = -s;
	rot[2][1] = s;
	return rot;
}

function RotateY(angle)
{
	var rot = matrix(3,3,0.0);
	var s = Math.sin(angle);
	var c = Math.cos(angle);
	rot[1][1] = 1.0;
	rot[2][2] = c;
	rot[0][0] = c;
	rot[2][0] = -s;
	rot[0][2] = s;
	return rot;
}
	
function RotateZ(angle)
{
	var rot = matrix(3,3,0.0);
	var s = Math.sin(angle);
	var c = Math.cos(angle);
	rot[2][2] = 1.0;
	rot[0][0] = c;
	rot[1][1] = c;
	rot[0][1] = -s;
	rot[1][0] = s;
	return rot;
}
	
function matrixMult(m, matrixArr)
{
	var result = matrix(3,3,0.0);
	result[0][0] = m[0][0] * matrixArr[0][0] + m[0][1] * matrixArr[1][0] + m[0][2] * matrixArr[2][0];
	result[0][1] = m[0][0] * matrixArr[0][1] + m[0][1] * matrixArr[1][1] + m[0][2] * matrixArr[2][1];
	result[0][2] = m[0][0] * matrixArr[0][2] + m[0][1] * matrixArr[1][2] + m[0][2] * matrixArr[2][2];
	result[1][0] = m[1][0] * matrixArr[0][0] + m[1][1] * matrixArr[1][0] + m[1][2] * matrixArr[2][0];
	result[1][1] = m[1][0] * matrixArr[0][1] + m[1][1] * matrixArr[1][1] + m[1][2] * matrixArr[2][1];
	result[1][2] = m[1][0] * matrixArr[0][2] + m[1][1] * matrixArr[1][2] + m[1][2] * matrixArr[2][2];
	result[2][0] = m[2][0] * matrixArr[0][0] + m[2][1] * matrixArr[1][0] + m[2][2] * matrixArr[2][0];
	result[2][1] = m[2][0] * matrixArr[0][1] + m[2][1] * matrixArr[1][1] + m[2][2] * matrixArr[2][1];
	result[2][2] = m[2][0] * matrixArr[0][2] + m[2][1] * matrixArr[1][2] + m[2][2] * matrixArr[2][2];
	return result;
}

function rotateVector(p, rot, x, y, z) {
    p[0] = rot[0][0]*x + rot[0][1]*y + rot[0][2]*z;
    p[1] = rot[1][0]*x + rot[1][1]*y + rot[1][2]*z;
    p[2] = rot[2][0]*x + rot[2][1]*y + rot[2][2]*z;
}