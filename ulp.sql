PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE IF NOT EXISTS "position_measurement" (
	"id"	INTEGER NOT NULL,
	"ra"	REAL NOT NULL,
	"dec"	REAL NOT NULL,
	"err_ellipse_r0"	REAL,
	"err_ellipse_r1"	REAL,
	"err_ellipse_ang"	INTEGER,
	"ulp_id"	INTEGER NOT NULL,
	FOREIGN KEY("ulp_id") REFERENCES "ulp"("id"),
	PRIMARY KEY("id" AUTOINCREMENT)
);
CREATE TABLE IF NOT EXISTS "user" (
	"id"	INTEGER NOT NULL,
	"username"	TEXT NOT NULL,
	"email"	TEXT NOT NULL UNIQUE,
	PRIMARY KEY("id" AUTOINCREMENT)
);
CREATE TABLE IF NOT EXISTS "collaboration" (
	"id"	INTEGER NOT NULL,
	"name"	TEXT NOT NULL UNIQUE,
	PRIMARY KEY("id" AUTOINCREMENT)
);
CREATE TABLE IF NOT EXISTS "user_collaboration" (
	"id"	INTEGER NOT NULL,
	"user_id"	INTEGER NOT NULL,
	"collaboration_id"	INTEGER NOT NULL,
	FOREIGN KEY("collaboration_id") REFERENCES "collaboration"("id"),
	FOREIGN KEY("user_id") REFERENCES "user"("id"),
	PRIMARY KEY("id" AUTOINCREMENT)
);
CREATE TABLE IF NOT EXISTS "ulp" (
	"id"	INTEGER NOT NULL,
	"name"	TEXT UNIQUE,
	"uploader_id"	INTEGER,
	"upload_date"	TEXT,
	PRIMARY KEY("id" AUTOINCREMENT),
	FOREIGN KEY("uploader_id") REFERENCES "user"("id")
);
CREATE TABLE IF NOT EXISTS "ulp_permission" (
	"id"	INTEGER NOT NULL,
	"ulp_id"	INTEGER NOT NULL,
	"collaboration_id"	INTEGER NOT NULL,
	FOREIGN KEY("collaboration_id") REFERENCES "collaboration"("id"),
	PRIMARY KEY("id" AUTOINCREMENT),
	FOREIGN KEY("ulp_id") REFERENCES "ulp"("id")
);
CREATE TABLE IF NOT EXISTS "parameter" (
	"id"	INTEGER NOT NULL,
	"name"	TEXT NOT NULL UNIQUE,
	"description"	TEXT,
	PRIMARY KEY("id" AUTOINCREMENT)
);
CREATE TABLE IF NOT EXISTS "parameter_measurement" (
	"id"	INTEGER NOT NULL,
	"ulp_id"	INTEGER NOT NULL,
	"parameter_id"	INTEGER NOT NULL,
	"value"	TEXT NOT NULL,
	"error"	TEXT,
	FOREIGN KEY("ulp_id") REFERENCES "ulp"("id"),
	FOREIGN KEY("parameter_id") REFERENCES "parameter"("id"),
	PRIMARY KEY("id" AUTOINCREMENT)
);
DELETE FROM sqlite_sequence;
INSERT INTO sqlite_sequence VALUES('ulp',0);
COMMIT;
